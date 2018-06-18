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

function step!(solver::ContinuationSolver)
    predictor_corrector_step!(solver.cache, solver.opts)
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

function step!(solver::ContinuationSolver, max_steps)
    cache = solver.cache
    for _ in 1:max_steps
        step!(solver)
        if ! isindomain(cache.u, cache.prob_cache)
            return true
        end
    end
    return false
end

function sweep!(solver;
                u0 = get_u0(solver.cache.prob_cache.prob),
                h = solver.opts.h0,
                past_points = [],
                direction = solver.opts.direction,
                max_steps = solver.opts.max_samples)
    opts = solver.opts
    cache = solver.cache

    cache.direction = direction
    cache.u = u0
    new_sweep!(solver.sol, direction)
    for u in past_points
        push_point!(solver.sol, u)
    end
    push_point!(solver.sol, u0)
    step!(solver, max_steps)
end

function solve!(solver::ContinuationSolver)
    opts = solver.opts
    cache = solver.cache
    H = residual!(cache.H, cache.u, cache.prob_cache)
    if ! isalmostzero(H, opts.rtol, opts.atol)
        error("Initial point is not a root: norm(H) = ", norm(H))
        # cache.u = nearest_root!(cache.u, opts.rtol, opts.atol)
    end

    sweep!(solver)
    sweep!(solver; direction = solver.opts.direction * -1)

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
            sweep!(solver;
                   u0 = u1,
                   past_points = [u0],
                   direction = direction,
                   h = h)
            append!(bifurcations, solver.sol.sweeps[end].simple_bifurcation)
        end
    end

    return solver
end
