module TestJacobian
include("preamble.jl")

using StaticArrays: SVector
using SymEngine: diff, subs, symbols

using Bifurcations.Continuations: residual!, residual_jacobian!, get_prob_cache
using Bifurcations.Codim2: DiffEqCodim2Problem, NormalizingAS
using Bifurcations.Examples: PredatorPrey

sjac(f, x) = [diff(f[i], x[j]) for i in 1:length(f), j in 1:length(x)]

nx = length(PredatorPrey.u0)
np = length(PredatorPrey.p)
f = PredatorPrey.f
pe = [symbols("p_$i") for i in 1:np]
xe = [symbols("x_$i") for i in 1:nx]
ve = [symbols("v_$i") for i in 1:nx]
all_symbols = vcat(pe, xe, ve)
fe = f(xe, pe, 0)
fe = Vector(fe)  # just for safety
dfdx = sjac(fe, xe)

ue = vcat(xe, ve, pe[1], pe[3])
He = vcat(fe,
          dfdx * ve,
          ve ⋅ ve - 1)
dHdu = sjac(He, ue)

prob = DiffEqCodim2Problem(
    PredatorPrey.ode,
    (@lens _[1]),
    (@lens _[3]),
    ([-Inf, -Inf], [Inf, Inf]);
    augmented_system = NormalizingAS(),
)

function test_predator_prey(
        xr, vr, p1, p3;
        prob_cache = get_prob_cache(prob))

    H_ = nothing
    J_ = nothing

    prob = prob_cache.prob
    pr = collect(Float64, prob.de_prob.p)
    ur = SVector(vcat(xr, vr, pr[1], pr[3])...) :: SVector{6, Float64}
    all_values = vcat(pr, xr, vr) :: Vector{Float64}
    evaluated = (e) -> Float64(subs(e, zip(all_symbols, all_values)...))

    H_desired = evaluated.(He)
    J_desired = evaluated.(dHdu)

    H1_actual = residual!(H_, ur, prob_cache)
    @test H1_actual ≈ H_desired

    H2_actual, J_actual = residual_jacobian!(H_, J_, ur, prob_cache)
    @test H2_actual ≈ H_desired
    @test J_actual ≈ J_desired
end

rng = MersenneTwister(0)
@testset begin
    for args in [
            ([0, 0], [0, 1], 0.6, 5),
            ([0, 0], [-1, 0], 0.384806, 7.79613),
            ]
        test_predator_prey(args...)
    end
    for _ in 1:100
        p1 = randn(rng)
        p3 = randn(rng)
        xr = randn(rng, nx)
        vr = normalize(randn(rng, nx))
        test_predator_prey(xr, vr, p1, p3)
    end
end

end  # module
