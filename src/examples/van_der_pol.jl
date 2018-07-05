module VanDerPol

using DiffEqBase: ODEProblem
using Parameters: @with_kw, @unpack
using Setfield: @lens

using ...Bifurcations: BifurcationProblem

@with_kw struct VanDerPolParam{D, A, Ω}
    d::D = 1.0  # 5.0
    a::A = 0.0  # 5.0
    ω::Ω = 2.466
end

function f(du, u, p, t)
    @unpack d, a, ω = p
    du[1] = u[2]
    du[2] = -u[1] - d * (u[1]^2 - 1) * u[2] + a * cos(ω * t)
    nothing
end

make_prob(
        p = VanDerPolParam();
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
