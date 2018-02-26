"""
The normal form of pitchfork bifurcations.

https://en.wikipedia.org/wiki/Pitchfork_bifurcation
"""
module Pitchfork

using DiffEqBase: ODEProblem
using Setfield: @lens

using ...Bifurcations: FixedPointBifurcationProblem

f(u, p, t) = p * u - u^3

u0 = 0.0
tspan = (0.0, 1.0)
p = -1.0
ode = ODEProblem(f, u0, tspan, p)

param_axis = @lens _
prob = FixedPointBifurcationProblem(ode, param_axis, (-1.0, 1.0))

end  # module
