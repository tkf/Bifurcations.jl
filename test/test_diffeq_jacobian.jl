module TestDiffeqJacobian
include("preamble.jl")

using Bifurcations.Continuations: get_prob_cache, get_u0, residual_jacobian
using Bifurcations.Examples.Calcium

@testset "ODEFunction (inplace)" begin
    prob0 = Calcium.make_prob()
    prob1 = Calcium.make_prob(ode=Calcium.make_ode_jac())
    cache0 = get_prob_cache(prob0)
    cache1 = get_prob_cache(prob1)

    @test prob0.p.de_prob.f.jac === nothing
    @test prob1.p.de_prob.f.jac === Calcium.jac

    u0 = get_u0(prob0)
    rng = MersenneTwister(0)
    for _ in 1:5
        u = typeof(u0)(randn(rng, eltype(u0), length(u0)))

        H0, J0 = residual_jacobian(u, cache0)
        H1, J1 = residual_jacobian(Array(u), cache1)

        @test H1 ≈ H0
        @test J1 ≈ J0
    end
end

end  # module
