module TestCodim2
using Base.Test

using Bifurcations
using Bifurcations: BifurcationProblem, special_points
using Bifurcations.Examples: Calcium
using Setfield: @lens
include("utils.jl")

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
    end
end

end  # module
