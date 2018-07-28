"""
Modified Morris-Lecar model from [Dhooge, Govaerts, Kuznetsov (2003)].

* [Dhooge, Govaerts, Kuznetsov (2003)].
  Numerical Continuation of Fold Bifurcations of Limit Cycles in MATCONT

[Dhooge, Govaerts, Kuznetsov (2003)]: https://doi.org/10.1007/3-540-44860-8_72
"""
module MorrisLecar

using DiffEqBase: ODEProblem
using Parameters: @with_kw, @unpack
using StaticArrays: SVector
using Setfield: @lens

using ...Bifurcations: BifurcationProblem
using ...Codim2: DiffEqCodim2Problem

@with_kw struct MorrisLecarParam{yType, zType}
    y::yType = 0.1104711840729923
    z::zType = 0.1
end

function f(u, p::MorrisLecarParam, t)
    @unpack y, z = p
    v, w = u
    m∞ = (1 + tanh((v + 0.01) / 0.15)) / 2
    w∞ = (1 + tanh((v - z) / 0.145)) / 2
    τ = cosh((v - 0.1) / 0.29)
    return SVector(
        y - 0.5(v + 0.5) - 2w * (v + 0.7) - m∞ * (v - 1),
        1.15(w∞ - w) * τ,
    )
end

function f(du, u, p::MorrisLecarParam, t)
    du .= f(SVector{2}(u), p, t)
    nothing
end


function make_prob(
        p = MorrisLecarParam();
        u0 = SVector(0.047222144637300754, 0.32564026295750237),
        tspan = (0.0, 30.0),
        ode = ODEProblem{!(u0 isa SVector)}(f, u0, tspan, p),
        param_axis = (@lens _.y),
        t_domain = (-0.05, 0.2),
        kwargs...)
    return BifurcationProblem(
        ode, param_axis, t_domain;
        kwargs...)
end

prob = make_prob()
ode = prob.p.de_prob
p = ode.p
u0 = ode.u0
param_axis = prob.p.param_axis

end  # module
