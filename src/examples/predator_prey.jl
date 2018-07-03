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

u0 = SVector(0.0, 0.0)
tspan = (0.0, 30.0)
p = (0.0, 3, 5, 3)
ode = ODEProblem{false}(f, u0, tspan, p)

param_axis = @lens _[1]
prob = BifurcationProblem(ode, param_axis, (0.0, 1.0),
                          phase_space = (SVector(-0.2, -Inf),  # u_min
                                         SVector(1.0, Inf)))   # u_max

end  # module
