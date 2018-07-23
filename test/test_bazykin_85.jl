module TestBazykin85
include("preamble_plots.jl")

using DiffEqBase: ODEProblem
using StaticArrays: SVector

using Bifurcations: Codim1, Codim2, resolved_points
using Bifurcations.BifurcationsBase: contkind, HopfCont
using Bifurcations.Codim2: NormalizingAS, BackReferencingAS
using Bifurcations.Examples: Bazykin85

# Bifurcation points calculated in some version of Bifurcations.jl
KNOWN_POINTS = Dict(
    # Bazykin85.Bazykin85Param => points
    Bazykin85.Bazykin85Param(
        ϵ = 0.01,
    ) => [
        (Codim2.PointTypes.cusp,
         [32.5936, 20.4743, 0.956835, 0.290633, 0.901233, 0.00356842]),
        (Codim2.PointTypes.bogdanov_takens,
         [42.688, 21.6306, 0.993354, 0.1151, 0.860704, 0.00605879]),
        (Codim2.PointTypes.bogdanov_takens,
         [6.26577, 3.60155, 0.932668, 0.360736, 0.453624, 0.175128]),
        (Codim2.PointTypes.cusp,
         [25.8231, 2.44209, 0.999458, 0.0329046, 0.0887673, 2.80236]),
    ],
)

@testset "smoke Bazykin85 codim-2 (ϵ=$ϵ, $(nameof(ASType)))" for
        ϵ in [0.01, 0.001],
        ASType in [NormalizingAS, BackReferencingAS]

    prob = Bazykin85.make_prob(
        Bazykin85.Bazykin85Param(
            ϵ = ϵ,
        ),
        phase_space = (SVector(0.0, 0.0),  # u_min
                       SVector(Inf, Inf)), # u_max
    )
    ode = prob.p.de_prob

    codim1_solver = init(prob)
    solve!(codim1_solver)

    if ϵ > 0.001
        @test length(special_points(codim1_solver)) == 4
    else
        @test_broken length(special_points(codim1_solver)) == 4
    end

    sn_point, = sort!(
        special_points(codim1_solver, Codim1.PointTypes.saddle_node),
        by=p->p.u0[1])
    hopf_point, = special_points(codim1_solver, Codim1.PointTypes.hopf)

    hopf_solver1 = nothing
    sn_solver1 = nothing
    for point in [sn_point, hopf_point]
        codim2_prob = BifurcationProblem(
            point,
            codim1_solver,
            (@lens _.δ),
            (0.0, 10.0),
            augmented_system = ASType(),
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

        @testset "resolve $p" for p in [Codim2.PointTypes.bogdanov_takens,
                                        Codim2.PointTypes.fold_hopf,
                                        Codim2.PointTypes.hopf_hopf]
            @test_nothrow resolved_points(codim2_solver, p)
        end

        if contkind(codim2_solver) isa HopfCont
            hopf_solver1 = codim2_solver
        else
            sn_solver1 = codim2_solver
        end
    end

    @testset "switch to saddle-node" begin
        @assert hopf_solver1 !== nothing
        # Find the right-most Bogdanov-Takens bifurcation:
        point = first(sort(
            special_points(hopf_solver1,
                           Codim2.PointTypes.bogdanov_takens);
            by = p -> p.u0[end - 1],  # α
            rev = true))
        sn_prob = BifurcationProblem(
            point, hopf_solver1;
            augmented_system = ASType(),
        )
        sn_solver2 = init(sn_prob)
        @test_nothrow solve!(sn_solver2)
    end

    @testset "switch to hopf" begin
        @assert sn_solver1 !== nothing
        # Find the right-most Bogdanov-Takens bifurcation:
        point = first(sort(
            special_points(sn_solver1,
                           Codim2.PointTypes.bogdanov_takens);
            by = p -> p.u0[end - 1],  # α
            rev = true))
        hopf_prob = BifurcationProblem(
            point, sn_solver1;
            augmented_system = ASType(),
        )
        hopf_solver2 = init(hopf_prob)
        @test_nothrow solve!(hopf_solver2)

        bpoints = special_points(hopf_solver2, Codim2.PointTypes.bautin)
        @test length(bpoints) == 1
        @test_nothrow resolved_points(hopf_solver2, Codim2.PointTypes.bautin)
    end

    if haskey(KNOWN_POINTS, ode.p)
        desired_points = KNOWN_POINTS[ode.p]
        @testset "known points" begin
            # TODO: do this for all solvers (hopf_solver1 etc.)
            uniquified = sort(
                uniquify_points(resolved_points(sn_solver1)),
                by = point -> point.u[end],
            )
            @test length(uniquified) == length(desired_points)
            for (actual, desired) in zip(uniquified, desired_points)
                type_desired, u_desired = desired
                if ASType == NormalizingAS
                    @test actual.u ≈ u_desired  rtol=1e-4
                else
                    # Compare DS state and parameters; ignore the eigenvector
                    # TODO: don't hard-code offsets:
                    @test actual.u[1:2] ≈ u_desired[1:2]  rtol=1e-4
                    @test actual.u[end-1:end] ≈ u_desired[end-1:end]  rtol=1e-4
                end
                @test actual.point_type == type_desired
            end
        end
    end
end

end  # module
