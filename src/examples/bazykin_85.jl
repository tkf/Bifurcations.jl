module Bazykin85

using DiffEqBase: ODEProblem
using Parameters: @with_kw, @unpack
using StaticArrays: SVector
using Setfield: @lens
import Setfield

using ...Bifurcations: BifurcationProblem
using ...Codim2: DiffEqCodim2Problem

@with_kw struct Bazykin85Param{A, E, G, D}
    α::A = 0.1
    ϵ::E = 0.01  # ≪ 1 ?
    γ::G = 1.0
    δ::D = 0.1
end


function f(x::SVector, p, t)
    @unpack α, ϵ, γ, δ = p
    return SVector(
             x[1] - x[1] * x[2] / (1 + α * x[1]) - ϵ * x[1]^2,
        -γ * x[2] + x[1] * x[2] / (1 + α * x[1]) - δ * x[2]^2,
    )
end

function f(du, u, p, t)
    du .= f(SVector{2}(u), p, t)
    nothing
end


make_prob(
        p = Bazykin85Param();
        u0 = SVector(1 / p.ϵ, 0.0),
        tspan = (0.0, 30.0),
        ode = ODEProblem{!(u0 isa SVector)}(f, u0, tspan, p),
        param_axis = (@lens _.α),
        t_domain = (0.01, 1.5),
        kwargs...) =
    BifurcationProblem(ode, param_axis, t_domain;
                       phase_space = (SVector(-0.1, -0.1),  # u_min
                                      SVector(Inf, Inf)),   # u_max
                       kwargs...)

prob = make_prob()
ode = prob.p.de_prob
u0 = ode.u0
param_axis = prob.p.param_axis
t_domain = prob.t_domain

function make_codim2_prob(
        p = Bazykin85Param();
        u0 = SVector(1 / p.ϵ, 0.0),
        tspan = (0.0, 30.0),
        ode = ODEProblem{!(u0 isa SVector)}(f, u0, tspan, p),
        param_axis1 = (@lens _.α),
        param_axis2 = (@lens _.δ),
        t_domain = ([0.01, 0.0], [1.5, 10.0]),
        kwargs...)
    return DiffEqCodim2Problem(
        ode,
        param_axis1,
        param_axis2,
        t_domain;
        kwargs...)
end

end  # module
