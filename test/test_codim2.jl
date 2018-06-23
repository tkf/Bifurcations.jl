module TestCodim2
using Base.Test

import GR  # Workaround Plots.jl world age problem
using Plots

using Bifurcations
using Bifurcations: plot  # TODO: stop doing this
using Bifurcations: BifurcationProblem, special_points
using Bifurcations.Examples: Calcium
using Setfield: @lens
include("utils.jl")

function nullshow(plt::Plots.Plot)
    nullshow(MIME("image/png"), plt)
end

@testset "smoke Calcium codim-2" begin
    codim1 = init(Calcium.prob)
    solve!(codim1)

    for point in special_points(codim1)
        codim2_prob = BifurcationProblem(
            point,
            codim1,
            (@lens _.gca),
            (0.0, 8.0),
        )
        codim2_solver = init(codim2_prob)
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
