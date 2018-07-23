"""
2D Predator-Prey Model

Taken from:
[14.3 pp2: A 2D Predator-Prey Model][auto07p/14.3]
of AUTO-07p manual.

[auto07p/14.3]: http://www.dam.brown.edu/people/sandsted/auto/auto07p.pdf#page=138
"""
module PredatorPrey

using DiffEqBase: ODEProblem
using StaticArrays: SVector
using Setfield: @lens
import Setfield

using ...Bifurcations: BifurcationProblem
using ...Codim2: cast_container, DiffEqCodim2Problem

f1(u, p) =
    p[2] * u[1] * (1 - u[1]) - u[1] * u[2] - p[1] * (1 - exp(- p[3] * u[1]))
f2(u, p) =
    - u[2] + p[4] * u[1] * u[2]

f(u::SVector, p, t) = SVector(f1(u, p), f2(u, p))

function f(du, u, p, t)
    du[1] = f1(u, p)
    du[2] = f2(u, p)
    nothing
end


function make_prob(
        p = (0.0, 3, 5, 3);
        u0 = SVector(0.0, 0.0),
        tspan = (0.0, 30.0),
        ode = ODEProblem{!(u0 isa SVector)}(f, u0, tspan, p),
        param_axis = (@lens _[1]),
        t_domain = (0.0, 1.0),
        kwargs...)
    x_min = SVector(-0.2, -Inf)
    x_max = SVector(1.0, Inf)
    x_min = cast_container(typeof(ode.u0), x_min)
    x_max = cast_container(typeof(ode.u0), x_max)
    return BifurcationProblem(
        ode, param_axis, t_domain;
        phase_space = (x_min, x_max),
        kwargs...)
end

prob = make_prob()
ode = prob.p.de_prob
p = ode.p
u0 = ode.u0
param_axis = prob.p.param_axis

function make_codim2_prob(
        p = ode.p;
        u0 = SVector(0.0, 0.0),
        tspan = (0.0, 30.0),
        ode = ODEProblem{!(u0 isa SVector)}(f, u0, tspan, p),
        param_axis1 = (@lens _[1]),
        param_axis2 = (@lens _[3]),
        t_domain = ([-Inf, -Inf], [Inf, Inf]),
        kwargs...)
    return DiffEqCodim2Problem(
        ode,
        param_axis1,
        param_axis2,
        t_domain;
        kwargs...)
end

end  # module
