module TestFoldLC
include("preamble.jl")

using Bifurcations: Codim1, Codim2, resolved_points, LimitCycleProblem,
    Codim1LimitCycle
using Bifurcations.Codim2LimitCycle: FoldLimitCycleProblem
using Bifurcations.Examples: Bautin

prob = Bautin.make_prob(Bautin.BautinParam(β₂ = 1.0))
solver0 = init(prob)
solve!(solver0)

hopf_points = resolved_points(solver0)
@test length(hopf_points) == 1
@test hopf_points[1].point_type === Codim1.PointTypes.hopf

prob_lc = LimitCycleProblem(
    hopf_points[1], solver0;
    num_mesh = 20,
    degree = 3,
)
solver_lc = init(
    prob_lc;
    start_from_nearest_root = true,
    bidirectional_first_sweep = false,
    max_branches = 0,
    nominal_angle_rad = 0.01,
)
solve!(solver_lc)

@test_broken length(solver_lc.sweeps[1].super.simple_bifurcation) == 0

flc_points = resolved_points(solver_lc)
@test_broken length(flc_points) == 1
@test flc_points[1].point_type === Codim1LimitCycle.PointTypes.saddle_node

#=
prob_flc = FoldLimitCycleProblem(
    flc_points[1], solver_lc;
)
=#

end  # module
