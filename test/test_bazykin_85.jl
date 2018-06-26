module TestBazykin85
include("preamble_plots.jl")

using Bifurcations.Examples: Bazykin85

@testset "smoke Bazykin85 codim-2" begin
    codim1_solver = init(Bazykin85.prob)
    solve!(codim1_solver)

    point_list = sort!(special_points(codim1_solver), by=p->p.u0[end])
    for point in point_list[2:3]  # TODO: avoid manual indexing
        codim2_prob = BifurcationProblem(
            point,
            codim1_solver,
            (@lens _.δ),
            (0.0, 10.0),
        )
        codim2_solver = init(
            codim2_prob;
            # h0 = 0.001,
            nominal_angle_rad = 0.01,
            max_samples = 1000,
        )
        solve!(codim2_solver)

        @testset "show" begin
            smoke_test_solver_show(codim2_solver)
        end

        @testset "plot" begin
            smoke_test_solver_plot(codim2_solver)
        end
    end
end


end  # module
