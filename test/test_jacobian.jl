module TestJacobian
include("preamble.jl")

using DiffEqBase: ODEProblem
using StaticArrays: SVector
using SymEngine: diff, subs, symbols

using Bifurcations.Continuations: residual, residual_jacobian, get_prob_cache
using Bifurcations.Codim2: DiffEqCodim2Problem, NormalizingAS
using Bifurcations.Examples: PredatorPrey

sjac(f, x) = [diff(f[i], x[j]) for i in 1:length(f), j in 1:length(x)]

nx = length(PredatorPrey.u0)
np = length(PredatorPrey.p)
f = generalized_f(PredatorPrey)
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


make_prob(;
        p = (0.0, 3, 5, 3),
        u0 = SVector(0.0, 0.0),
        tspan = (0.0, 30.0),
        ode = ODEProblem{!(u0 isa SVector)}(PredatorPrey.f, u0, tspan, p),
        kwargs...) =
    DiffEqCodim2Problem(
        ode,
        (@lens _[1]),
        (@lens _[3]),
        ([-Inf, -Inf], [Inf, Inf]);
        augmented_system = NormalizingAS(),
        kwargs...)

function test_predator_prey(
        xr, vr, p1, p3;
        u0 = SVector(0.0, 0.0),
        prob = make_prob(u0 = u0),
        prob_cache = get_prob_cache(prob))

    pr = collect(Float64, prob.de_prob.p)
    pr[1] = p1
    pr[3] = p3

    ur = vcat(xr, vr, pr[1], pr[3])
    if prob.de_prob.u0 isa SVector
        ur = SVector(ur...)
    else
        utype = typeof(prob.de_prob.u0)
        if !(ur isa utype)
            ur = utype(ur)
        end
    end

    all_values = vcat(pr, xr, vr) :: Vector{Float64}
    evaluated = (e) -> Float64(subs(e, zip(all_symbols, all_values)...))

    H_desired = evaluated.(He)
    J_desired = evaluated.(dHdu)

    H1_actual = residual(ur, prob_cache)
    @test H1_actual ≈ H_desired

    H2_actual, J_actual = residual_jacobian(ur, prob_cache)
    @test H2_actual ≈ H_desired
    @test J_actual ≈ J_desired
end

rng = MersenneTwister(0)
@testset begin
    for args in [
            ([0, 0], [0, 1], 0.6, 5),
            ([0, 0], [-1, 0], 0.384806, 7.79613),
            ],
        u0 in [SVector(0.0, 0.0), [0.0, 0.0]]

        test_predator_prey(args...; u0=u0)
    end
    for _ in 1:100,
        u0 in [SVector(0.0, 0.0), [0.0, 0.0]]

        p1 = randn(rng)
        p3 = randn(rng)
        xr = randn(rng, nx)
        vr = normalize(randn(rng, nx))
        test_predator_prey(xr, vr, p1, p3; u0=u0)
    end
end

end  # module
