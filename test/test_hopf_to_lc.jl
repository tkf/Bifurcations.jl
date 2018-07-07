module TestHopfToLC
include("preamble_plots.jl")

using Bifurcations
using Bifurcations: Codim1, special_points, LimitCycleProblem, limitcycles
using Bifurcations.Examples: PredatorPrey

solver0 = init(PredatorPrey.prob)
solve!(solver0)

hopf_point, = special_points(solver0, Codim1.PointTypes.hopf)

period_max = 36.0
prob_lc = LimitCycleProblem(
    hopf_point, solver0;
    period_bound = (0.0, period_max),
    num_mesh = 20,
    degree = 3,
)
solver1 = init(
    prob_lc;
    start_from_nearest_root = true,
)
solve!(solver1)

periods = [lc.period for lc in limitcycles(solver1)]
@test length(periods) > 0
@test all(!isnan, periods)
@test maximum(periods) >= period_max  # it is stopped due to out-of-domain
@test minimum(periods) > 0

plt_multi_panel = plot(
    plot(solver1, (:x=>:parameter, :y=>1, :color=>:period)),
    plot(solver1, (:x=>:parameter, :y=>:period, :color=>:period)),
    plot_state_space(solver1);
    clims = (10, period_max),  # TODO: make it work without this
)

@test_nothrow nullshow(plt_multi_panel)

end  # module
