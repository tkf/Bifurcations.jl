module TestPredatorPrey
using Base.Test

import GR  # Workaround Plots.jl world age problem
using Plots

using Bifurcations
using Bifurcations: plot  # TODO: stop doing this
using Bifurcations: BifurcationProblem, special_points
using Bifurcations.Examples: PredatorPrey
using Setfield: @lens
include("utils.jl")

function nullshow(plt::Plots.Plot)
    nullshow(MIME("image/png"), plt)
end

@testset "smoke PredatorPrey codim-2" begin
    codim1_solver = init(PredatorPrey.prob)
    solve!(codim1_solver)

    for point in special_points(codim1_solver)
        codim2_prob = BifurcationProblem(
            point,
            codim1_solver,
            (@lens _[3]),
            (0.0, 10.0),
        )
        # TODO: Implement detection of Bogdanov-Takens bifurcation and
        # stop manually setting max_branches=0.
        codim2_solver = init(codim2_prob; max_branches=0)
        @test_nothrow solve!(codim2_solver)

        @testset "show" begin
            smoke_test_solver_show(codim2_solver)
        end

        @testset "plot" begin
            smoke_test_solver_plot(codim2_solver)
        end
    end
end

end  # module
