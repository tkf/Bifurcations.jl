module TestCodim2
include("preamble_plots.jl")

using Bifurcations.Codim2: NormalizingAS, BackReferencingAS
using Bifurcations.Examples: Calcium

@testset "smoke Calcium codim-2 ($(nameof(ASType)))" for
        ASType in [NormalizingAS, BackReferencingAS]

    codim1 = init(Calcium.prob)
    solve!(codim1)

    for point in special_points(codim1)
        codim2_prob = BifurcationProblem(
            point,
            codim1,
            (@lens _.gca),
            (0.0, 8.0),
            augmented_system = ASType(),
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
