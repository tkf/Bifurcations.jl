"""
The normal form of transcritical bifurcations.

https://en.wikipedia.org/wiki/Transcritical_bifurcation
"""
module Transcritical

using DiffEqBase: ODEProblem
using Setfield: @lens

using ...Bifurcations: FixedPointBifurcationProblem

f(u, p, t) = p * u - u^2

u0 = 0.0
tspan = (0.0, 1.0)
p = -1.0
ode = ODEProblem(f, u0, tspan, p)

param_axis = @lens _
prob = FixedPointBifurcationProblem(ode, param_axis, (-1.0, 1.0))


using StaticArrays: @SVector

function closest_analytic(u, p)
    analytic = @SVector [0, p]
    _, i = findmin(abs.(u .- analytic))
    return [analytic[i], p]
end

end  # module
