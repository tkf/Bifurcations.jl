module Cusp

using DiffEqBase: ODEProblem
using Parameters: @with_kw, @unpack
using StaticArrays: SVector
using Setfield: @lens

using ...Bifurcations: BifurcationProblem

@with_kw struct CuspParam{B1, B2, S}
    # parameters are chosen so that u=1 is a fixed point.
    β₁::B1 = 1 / 2
    β₂::B2 = 1 / 2
    s::S = -1
end

function f(u::Real, p, t)
    @unpack β₁, β₂, s = p
    return β₁ + β₂ * u + s * u^3
end
# http://www.scholarpedia.org/article/Cusp_bifurcation#One-dimensional_Case

f(u::SVector, p, t) = SVector(f(u[1], p, t))


make_prob(
        p = CuspParam();
        u0 = SVector(1.0),
        tspan = (0.0, 30.0),
        ode = ODEProblem(f, u0, tspan, p),
        param_axis = (@lens _.β₁),
        t_domain = (-2.0, 2.0),
        kwargs...) =
    BifurcationProblem(ode, param_axis, t_domain;
                       kwargs...)

prob = make_prob()
ode = prob.p.de_prob
u0 = ode.u0
param_axis = prob.p.param_axis

end  # module
