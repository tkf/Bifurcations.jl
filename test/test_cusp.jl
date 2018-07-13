module TestCusp
include("preamble_plots.jl")

using Bifurcations: Codim1, Codim2, resolved_points

solver1 = init(Bifurcations.Examples.Cusp.prob)
solve!(solver1)

codim1_points = resolved_points(solver1)
@test length(codim1_points) == 2
@test codim1_points[1].point_type === Codim1.PointTypes.saddle_node
@test codim1_points[2].point_type === Codim1.PointTypes.saddle_node

sn_prob = BifurcationProblem(
    codim1_points[1],
    solver1,
    (@lens _.β₂),
    (-2.0, 2.0),
)
sn_solver = init(
    sn_prob;
)
solve!(sn_solver)

β₁s = [u[end-1] for sweep in sn_solver.super.sol.sweeps for u in sweep.u]
β₂s = [u[end]   for sweep in sn_solver.super.sol.sweeps for u in sweep.u]
@test all(@. isapprox(β₁s^2, (2 / (3 * √3) *  β₂s^(3 / 2))^2; atol=1e-6))
@test maximum(β₁s) > 1e-6
@test minimum(β₁s) < -1e-6
@test maximum(β₂s) > 2
@test minimum(β₂s) > -1e-6

codim2_points = resolved_points(sn_solver)
@test length(codim2_points) == 1
@test codim2_points[1].point_type === Codim2.PointTypes.cusp
β_cusp = codim2_points[1].u[end-1:end]
@test all(@. abs(β_cusp) < 1e-6)

plt = plot(sn_solver)
@test_nothrow nullshow(plt)

end  # module
