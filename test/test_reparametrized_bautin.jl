module TestReparametrizedBautin
include("preamble.jl")

using ForwardDiff
using Setfield: compose
using StaticArrays: SVector, SMatrix

using Bifurcations: Codim1, Codim2, resolved_points, reparametrize
using Bifurcations.Codim2: first_lyapunov_coefficient
using Bifurcations.Codim2LimitCycle: FoldLimitCycleProblem
using Bifurcations.Continuations: get_prob_cache, get_u0,
    residual!, residual_jacobian
using Bifurcations.Examples: Bautin
using Bifurcations.Examples.Reparametrization: orig_p

# Putting @info so that Travis won't stop testing.
@testset for opts in [
        [],
        [:extra => SVector(-1.0)],
        [:extra => SVector(-1.0), :shift => 0.1, :seed => 1],
        ]
    @info "TestReparametrizedBautin:"
    @info "opts = $opts"

    prob = reparametrize(
        Bautin.make_prob(),
        ; opts...)
    solver1 = init(prob)
    solve!(solver1)

    codim1_points = resolved_points(solver1)
    @test length(codim1_points) == 1
    @test codim1_points[1].point_type === Codim1.PointTypes.hopf
    @test codim1_points[1].u[end] ≈ 0  atol=1e-6

    @info "Continuation of Hopf bifurcation..."
    hopf_prob = BifurcationProblem(
        codim1_points[1],
        solver1,
        compose((@lens _.β₂), orig_p),
        (-2.0, 2.0),
    )
    hopf_solver = init(
        hopf_prob;
        verbose = true,
    )
    @time solve!(hopf_solver)

    hopf_β₁ = [u[end-1] for sweep in hopf_solver.super.sol.sweeps for u in sweep.u]
    hopf_β₂ = [u[end]   for sweep in hopf_solver.super.sol.sweeps for u in sweep.u]
    @test all(@. abs(hopf_β₁) < 1e-6)
    @test maximum(hopf_β₂) > 2
    @test minimum(hopf_β₂) < -2

    l1s = [
        first_lyapunov_coefficient(hopf_solver.cache.super.prob_cache, u)
        for sweep in hopf_solver.super.sol.sweeps for u in sweep.u
    ]
    if !all(@. (l1s > 0) == (hopf_β₂ > 0))
        println("[hopf_β₂ l1s]' =")
        display([hopf_β₂ l1s]')
    end
    @test all(@. (l1s > 0) == (hopf_β₂ > 0))

    @info "Resolving Bautin point..."
    @time codim2_points = resolved_points(hopf_solver)
    @test length(codim2_points) == 1
    @test codim2_points[1].point_type === Codim2.PointTypes.bautin
    β_bautin = codim2_points[1].u[end-1:end]
    @test all(@. abs(β_bautin) < 1e-6)

    @info "Creating FoldLimitCycleProblem..."
    @time flc_prob = FoldLimitCycleProblem(
        codim2_points[1],
        hopf_solver;
        num_mesh = 20,
        degree = 3,
    )
    @show flc_prob.super.num_mesh
    @show flc_prob.super.degree
    @test flc_prob.t_domain == ([-2.0, -2.0], [2.0, 2.0])

    @info "Checking Jacobian..."
    @testset "Manual Jacobian vs ForwardDiff" begin
        prob_cache = get_prob_cache(flc_prob)
        u0 = get_u0(flc_prob)
        rng = MersenneTwister(0)
        for _ in 1:5
            # Trying to avoid state values whose abs is larger than 1;
            # fitting polynomial to random numbers is not nice thing
            # to try here...
            u_lc = 0.01 .* tanh.(randn(rng, eltype(u0), length(u0)))
            d = length(flc_prob.xs0)
            u_lc[d + 1] = abs(u_lc[d + 1])  # period > 0
            H_actual, J_actual = residual_jacobian(u_lc, prob_cache)

            H_desired = similar(@view u_lc[1:end-1])
            J_desired = similar(H_desired, (length(H_desired), length(u_lc)))
            ForwardDiff.jacobian!(
                J_desired,
                (y, x) -> residual!(y, x, prob_cache),
                H_desired,
                u_lc,
            )

            @test J_actual ≈ J_desired
            @test H_actual ≈ H_desired
        end
    end

    @info "Continuation of fold bifurcation of a limit cycle..."
    flc_solver = init(
        flc_prob;
        start_from_nearest_root = true,
        max_branches = 0,  # TODO: stop manually doing this
        verbose = true,
    )
    @time solve!(flc_solver)

    flc_β₁ = [u[end-1] for sweep in flc_solver.super.sol.sweeps for u in sweep.u]
    flc_β₂ = [u[end]   for sweep in flc_solver.super.sol.sweeps for u in sweep.u]
    @show maximum(@. abs(4 * flc_β₁ + flc_β₂^2))
    @test all(@. abs(4 * flc_β₁ + flc_β₂^2) < 5e-2)
    @test maximum(flc_β₂) > 2
    @test minimum(flc_β₂) > -5e-2
end

end  # module
