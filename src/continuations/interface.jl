using Parameters: @with_kw

function get_prob_cache end
function get_u0 end
function residual! end
function residual_jacobian! end
function isindomain end

abstract type AbstractContinuationProblem{iip} end
abstract type AbstractProblemCache{P} end

mutable struct ContinuationCache{PC <: AbstractProblemCache,
                                 uType, HType, JType, hType}
    prob_cache::PC
    u::uType
    H::HType
    J::JType
    h::hType
    direction::Int
    corrector_success::Bool
    adaptation_success::Bool
    simple_bifurcation::Bool
end

function ContinuationCache(prob_cache::AbstractProblemCache,
                           h::Real, direction::Int = 1)
    u0 = get_u0(prob_cache.prob)
    N = length(u0)
    return ContinuationCache(
        prob_cache,
        u0,
        _similar(u0, N - 1),     # H
        _similar(u0, N - 1, N),  # J
        h,
        direction,
        false,
        false,
        false,
    )
end

ContinuationCache(prob::AbstractContinuationProblem, args...) =
    ContinuationCache(get_prob_cache(prob), args...)


@with_kw struct ContinuationOptions
    direction::Int = 1
    h0::Float64 = 1.0
    h_min::Float64 = 1e-6
    rtol::Float64 = 0.01
    atol::Float64 = 1e-6
    max_samples::Int = 100
    max_adaptations::Int = 100
    max_corrector_steps::Int = 100
    nominal_contraction::Float64 = 0.8
    nominal_distance::Float64 = 0.1
    nominal_angle_rad::Float64 = 2π * (30 / 360)
end


mutable struct ContinuationSolver{P <: AbstractContinuationProblem,
                                  C <: ContinuationCache,
                                  S <: ContinuationSolution}
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

function solve!(solver::ContinuationSolver)
    opts = solver.opts
    cache = solver.cache
    H = residual!(cache.H, cache.u, cache.prob_cache)
    if ! isalmostzero(H, opts.rtol, opts.atol)
        error("Initial point is not a root: norm(H) = ", norm(H))
        # cache.u = nearest_root!(cache.u, opts.rtol, opts.atol)
    end

    step!(solver, solver.opts.max_samples)

    # Flip the direction and solve it in the opposite direction:
    cache.direction *= -1
    cache.u = get_u0(cache.prob_cache.prob)
    new_sweep!(solver.sol)
    push_point!(solver.sol, cache.u)
    step!(solver, solver.opts.max_samples)

    # TODO: Detect the case that the solution is isomorphic to the
    # circle.

    return solver
end
