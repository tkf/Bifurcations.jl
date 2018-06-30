module TestBazykin85
include("preamble_plots.jl")

using DiffEqBase: ODEProblem
using StaticArrays: SVector

using Bifurcations: Codim1, Codim2, resolved_points
using Bifurcations.BifurcationsBase: contkind, HopfCont
using Bifurcations.Examples: Bazykin85

@testset "smoke Bazykin85 codim-2 (ϵ=$ϵ)" for ϵ in [0.01, 0.001]

    ode = let
        p = Bazykin85.Bazykin85Param(
            ϵ = ϵ,
        )
        u0 = SVector(1 / p.ϵ, 0.0)
        ODEProblem(Bazykin85.f, u0, Bazykin85.tspan, p)
    end
    prob = BifurcationProblem(
        ode,
        Bazykin85.param_axis,
        (0.01, 1.5);
        phase_space = (SVector(0.0, 0.0),  # u_min
                       SVector(Inf, Inf)), # u_max
    )

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
        sn_prob = BifurcationProblem(point, hopf_solver1)
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
        hopf_prob = BifurcationProblem(point, sn_solver1)
        hopf_solver2 = init(hopf_prob)
        @test_nothrow solve!(hopf_solver2)
    end
end

end  # module
