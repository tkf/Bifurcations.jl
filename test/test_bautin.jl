module TestBautin
include("preamble_plots.jl")

using Bifurcations: Codim1, Codim2, resolved_points

solver1 = init(Bifurcations.Examples.Bautin.prob)
solve!(solver1)

codim1_points = resolved_points(solver1)
@test length(codim1_points) == 1
@test codim1_points[1].point_type === Codim1.PointTypes.hopf

hopf_prob = BifurcationProblem(
    codim1_points[1],
    solver1,
    (@lens _.β₂),
    (-2.0, 2.0),
)
hopf_solver = init(
    hopf_prob;
)
solve!(hopf_solver)

β₁s = [u[end-1] for sweep in hopf_solver.super.sol.sweeps for u in sweep.u]
β₂s = [u[end]   for sweep in hopf_solver.super.sol.sweeps for u in sweep.u]
@test all(@. abs(β₁s) < 1e-6)
@test maximum(β₂s) > 2
@test minimum(β₂s) < -2

codim2_points = resolved_points(hopf_solver)
@test length(codim2_points) == 1
@test codim2_points[1].point_type === Codim2.PointTypes.bautin

plt = plot(hopf_solver)
@test_nothrow nullshow(plt)

end  # module
