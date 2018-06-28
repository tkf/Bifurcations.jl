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
