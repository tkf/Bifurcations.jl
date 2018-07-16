using Parameters: @unpack

"""
    as(self, ::Type{T}) :: T

Manual Go-style type [embedding] as a replacement of "inheritance from
concrete type".

[embedding]: https://golang.org/doc/effective_go.html?#embedding
"""
as(self, ::Type{T}) where {T} = as(self.super, T) :: T
as(self::T, ::Type{T}) where {T} = self


abstract type AbstractContinuationSolver end

mutable struct ContinuationSolver{P <: AbstractContinuationProblem,
                                  C <: ContinuationCache,
                                  S <: ContinuationSolution,
                                  } <: AbstractContinuationSolver
    prob::P
    opts::ContinuationOptions
    cache::C
    sol::S
    i::Int
end

function ContinuationSolver(prob, opts)
    cache = ContinuationCache(prob, opts.h0, opts.direction)
    sol = ContinuationSolution(typeof(cache.u))
    return ContinuationSolver(prob, opts, cache, sol, 0)
end

init(prob::AbstractContinuationProblem; kwargs...) =
    ContinuationSolver(prob, ContinuationOptions(; kwargs...))

solve(prob::AbstractContinuationProblem; kwargs...) =
    solve!(init(prob; kwargs...)).sol

function step!(solver::ContinuationSolver)
    predictor_corrector_step!(solver.cache, solver.opts)

    # Errors are now thrown in predictor_corrector_step!.  I need to
    # reconsider the error handling...
    if ! solver.cache.adaptation_success
        error("Failed to adapt steplength h.")
    end
    if ! solver.cache.corrector_success
        error("Failed in corrector loop.")
    end

    record!(solver.sol, solver.cache)
    solver.i += 1
end

function record!(sol, cache)
    push_point!(sol, cache)
end

function step!(wrapper::AbstractContinuationSolver, max_steps)
    solver = as(wrapper, ContinuationSolver)
    cache = solver.cache
    cache.h = solver.opts.h0
    for _ in 1:max_steps
        step!(wrapper)
        if ! isindomain(cache.u, cache.prob_cache)
            return true
        end
    end
    return false
end

struct SweepSetup{uType}
    direction::Int
    u0::uType
    past_points::Array{uType}
    max_steps::Int
end

SweepSetup(solver::ContinuationSolver;
           direction = solver.opts.direction,
           u0 = get_u0(solver.cache.prob_cache.prob),
           past_points = [],
           max_steps = solver.opts.max_samples,
           h = nothing,  # currently ignored  # TODO: use it?
           ) =
    SweepSetup(direction, u0, past_points, max_steps)

SweepSetup(solver::AbstractContinuationSolver; kwargs...) =
    SweepSetup(as(solver, ContinuationSolver); kwargs...)

function new_sweep!(solver::ContinuationSolver, setup::SweepSetup)
    @unpack direction, u0, past_points = setup
    cache = solver.cache
    cache.direction = direction
    cache.u = u0

    new_sweep!(solver.sol, direction)
    for u in past_points
        push_point!(solver.sol, u)
    end
    push_point!(solver.sol, u0)
end

function sweep!(solver::AbstractContinuationSolver; kwargs...)
    setup = SweepSetup(solver; kwargs...)
    new_sweep!(solver, setup)
    step!(solver, setup.max_steps)
end

function solve!(wrapper::AbstractContinuationSolver)
    solver = as(wrapper, ContinuationSolver)
    opts = solver.opts
    cache = solver.cache
    H = residual!(cache.H, cache.u, cache.prob_cache)
    if ! isalmostzero(H, opts.rtol, opts.atol)
        if opts.start_from_nearest_root
            cache.u = nearest_root!(cache, opts)
        else
            error("Initial point is not a root: norm(H) = ", norm(H))
        end
    end

    u0 = copy(cache.u)
    sweep!(wrapper; u0=u0)
    sweep!(wrapper; u0=u0, direction = solver.opts.direction * -1)

    # TODO: Detect the case that the solution is isomorphic to the
    # circle.

    bifurcations = vcat(solver.sol.sweeps[end-1].simple_bifurcation,
                        solver.sol.sweeps[end].simple_bifurcation)
    for _ in 1:opts.max_branches
        if isempty(bifurcations)
            break
        end
        sbint = shift!(bifurcations)
        for (u0, u1, direction, h) in new_branches!(cache, opts, sbint)
            if ! (isindomain(u0, cache.prob_cache) &&
                  isindomain(u1, cache.prob_cache))
                # Stepped outside the domain.  Skip it.
                continue
            end
            sweep!(wrapper;
                   u0 = u1,
                   past_points = [u0],
                   direction = direction,
                   h = h)
            append!(bifurcations, solver.sol.sweeps[end].simple_bifurcation)
        end
    end

    return wrapper
end


function residual(u::AbstractArray, cache::AbstractProblemCache)
    H = similar(@view u[1:end-1])
    return residual!(H, u, cache)
end

function residual(u::SVector, cache::AbstractProblemCache)
    return residual!(nothing, u, cache)
end


function residual_jacobian(u::AbstractArray, cache::AbstractProblemCache)
    H = similar(@view u[1:end-1])
    J = similar(H, (length(H), length(u)))
    return residual_jacobian!(H, J, u, cache)
end

function residual_jacobian(u::SVector, cache::AbstractProblemCache)
    return residual_jacobian!(nothing, nothing, u, cache)
end


function residual_jacobian!(cache::ContinuationCache, u::AbstractArray)
    @unpack H, J, prob_cache = cache
    H, J = residual_jacobian!(H, J, u, prob_cache)
    cache.u = u
    cache.H = H
    cache.J = J
    return H, J
end
