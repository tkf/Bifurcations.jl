special_points(solver::BifurcationSolver, args...) =
    special_points(solver.sol, args...)

special_points(sol::BifurcationSolution, args...) =
    [point for sweep in sol.sweeps for point in special_points(sweep, args...)]

special_points(sweep::BifurcationSweep) = sweep.special_points
special_points(sweep::BifurcationSweep, point_type::Enum) =
    special_points(sweep, (point_type,))
special_points(sweep::BifurcationSweep,
               point_types::NTuple{N, <:Enum}) where N =
    filter(p -> p.point_type in point_types, special_points(sweep))


"""
    problem_of(ctx) :: BifurcationProblem

Get the problem associated with `ctx` (a solver, solution, or sweep).
"""
problem_of(solver::ContinuationSolver) = solver.prob
problem_of(sol::ContinuationSolution) = sol.prob
problem_of(sweep::ContinuationSweep) = problem_of(sweep.sol.value)

problem_of(solver::BifurcationSolver) = solver.prob
problem_of(sol::BifurcationSolution) = problem_of(as(sol, ContinuationSolution))
problem_of(sweep::BifurcationSweep) = problem_of(as(sweep, ContinuationSweep))
