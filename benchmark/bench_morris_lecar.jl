using BenchmarkTools

using Setfield: @lens

using Bifurcations
using Bifurcations: Codim1, Codim2, special_intervals
using Bifurcations.Continuations: as, pre_solve!, sweep!
using Bifurcations.Codim2LimitCycle: FoldLimitCycleProblem
using Bifurcations.Examples: MorrisLecar

function make_morris_lecar_hopf_solver()
    solver1 = init(MorrisLecar.make_prob())
    solve!(solver1)
    hopf_point, = special_intervals(solver1, Codim1.PointTypes.hopf)

    hopf_prob = BifurcationProblem(
        hopf_point,
        solver1,
        (@lens _.z),
        (-1.0, 1.0),
    )
    return init(
        hopf_prob;
        nominal_angle_rad = 0.01,
    )
end

function make_morris_lecar_flc_solver()
    hopf_solver = make_morris_lecar_hopf_solver()
    solve!(hopf_solver)
    bautin_point, = special_intervals(hopf_solver, Codim2.PointTypes.bautin)

    flc_prob = FoldLimitCycleProblem(
        bautin_point,
        hopf_solver;
        num_mesh = 30,
        degree = 4,
    )
    return init(
        flc_prob;
        start_from_nearest_root = true,
        bidirectional_first_sweep = false,
        max_branches = 0,
        max_samples = 3,
        nominal_angle_rad = 2π * (5 / 360),
    )
end

suite = BenchmarkGroup()
suite["sweep!(hopf_solver)"] = @benchmarkable(
    sweep!(hopf_solver),
    setup = (hopf_solver = pre_solve!(make_morris_lecar_hopf_solver())))
suite["sweep!(flc_solver)"] = @benchmarkable(
    sweep!(flc_solver),
    setup = (flc_solver = pre_solve!(make_morris_lecar_flc_solver())))

if !isdefined(:SUITE)
    SUITE = BenchmarkGroup()
end
SUITE["morris_lecar"] = suite
