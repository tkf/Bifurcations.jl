module BazykinKhibnik81

using DiffEqBase: ODEProblem
using Parameters: @with_kw, @unpack
using StaticArrays: SVector
using Setfield: @lens
import Setfield

using ...Bifurcations: BifurcationProblem
using ...Codim2: DiffEqCodim2Problem

@with_kw struct BazykinKhibnik81Param{N, M, G}
    n::N = 0.5
    m::M = 0.5
    γ::G = 1
end

function f(u::SVector, p::BazykinKhibnik81Param, t)
    @unpack n, m, γ = p
    x, y = u
    return SVector(
        x * x * (1 - x) / (n + x) - x * y,
        - γ * y * (m - x),
    )
end

function f(du, u, p, t)
    du .= f(SVector{2}(u), p, t)
    nothing
end


make_prob(
        p = BazykinKhibnik81Param();
        u0 = SVector(1.0, 0.0),
        tspan = (0.0, 30.0),
        ode = ODEProblem{!(u0 isa SVector)}(f, u0, tspan, p),
        param_axis = (@lens _.m),
        t_domain = (-1e-3, 2.0),
        kwargs...) =
    BifurcationProblem(ode, param_axis, t_domain;
                       kwargs...)

prob = make_prob()
ode = prob.p.de_prob
u0 = ode.u0
param_axis = prob.p.param_axis
t_domain = prob.t_domain

function make_codim2_prob(
        p = BazykinKhibnik81Param();
        u0 = SVector(1.0, 0.0),
        tspan = (0.0, 30.0),
        ode = ODEProblem{!(u0 isa SVector)}(f, u0, tspan, p),
        param_axis1 = (@lens _.m),
        param_axis2 = (@lens _.n),
        t_domain = ([-1e-3, -1e-3], [2.0, 2.0]),
        kwargs...)
    return DiffEqCodim2Problem(
        ode,
        param_axis1,
        param_axis2,
        t_domain;
        kwargs...)
end

end  # module
