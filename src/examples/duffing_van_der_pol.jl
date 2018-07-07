module DuffingVanDerPol

using DiffEqBase: ODEProblem
using Parameters: @with_kw, @unpack
using Setfield: @lens

using ...Bifurcations: BifurcationProblem

@with_kw struct DuffingVanDerPolParam{A, B, D}
    a::A = 1.0
    b::B = 0.0
    d::D = 1.0
end

function f(du, u, p, t)
    @unpack a, b, d = p
    du[1] = u[2]
    du[2] = - a * u[1] - b * u[1]^3 - d * (u[1]^2 - 1) * u[2]
    nothing
end

make_prob(
        p = DuffingVanDerPolParam();
        u0 = [0.0, 0.0],
        tspan = (0.0, 2π),
        ode = ODEProblem(f, u0, tspan, p),
        param_axis = (@lens _.d),
        t_domain = (-2.0, 8.0),
        kwargs...) =
    BifurcationProblem(ode, param_axis, t_domain;
                       kwargs...)

prob = make_prob()
ode = prob.p.de_prob
param_axis = prob.p.param_axis
t_domain = prob.t_domain

end  # module
