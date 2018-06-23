special_points(solver::BifurcationSolver) = special_points(solver.sol)
special_points(sol::BifurcationSolution) =
    [point for sweep in sol.sweeps for point in special_points(sweep)]
special_points(sweep::BifurcationSweep) = sweep.special_points
