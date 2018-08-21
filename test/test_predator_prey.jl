module TestPredatorPrey
include("preamble_plots.jl")

using Bifurcations: resolved_points
using Bifurcations.Codim1
using Bifurcations.Codim2: NormalizingAS, BackReferencingAS
using Bifurcations.Examples.PredatorPrey

# Bifurcation points calculated in some version of Bifurcations.jl
CODIM1_KNOWN_POINTS = Dict(
    # PredatorPrey parameter (tuple) => points
    (0.0, 3, 5, 3) => [
        (Codim1.PointTypes.saddle_node,
         [0.0, 0.0, 0.6]),
        (Codim1.PointTypes.saddle_node,
         [1/3, 0, 0.821905296656178]),
        (Codim1.PointTypes.saddle_node,
         [0.41116155950979133, 0, 0.8329293250960671]),
        (Codim1.PointTypes.hopf,
         [1/3, 0.36576153597718997, 0.6715938474788321]),
    ],
)

@testset "smoke PredatorPrey" begin

    codim1_solver = init(PredatorPrey.prob, nominal_angle_rad=0.01)
    solve!(codim1_solver)

    @testset "known codim1 points" begin
        resolved = resolved_points(codim1_solver)
        known = CODIM1_KNOWN_POINTS[codim1_solver.prob.p.de_prob.p]
        @test length(resolved) == length(known)
        for (actual, desired) in zip(resolved, known)
            type_desired, u_desired = desired
            @test actual.point_type == type_desired
            @test actual.u ≈ u_desired  rtol=1e-4
        end
    end

    @testset "codim-2 ($(nameof(ASType))) from $(point.u0[[1, end]])" for
            ASType in [NormalizingAS, BackReferencingAS],
            point in special_intervals(codim1_solver)

        codim2_prob = BifurcationProblem(
            point,
            codim1_solver,
            (@lens _[3]),
            (0.0, 10.0),
            augmented_system = ASType(),
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
