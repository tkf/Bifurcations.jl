"""
The normal form of pitchfork bifurcations.

https://en.wikipedia.org/wiki/Pitchfork_bifurcation
"""
module Pitchfork

using DiffEqBase: ODEProblem
using Setfield: @lens

using ...Bifurcations: BifurcationProblem

f(u, p, t) = p * u - u^3

u0 = 0.0
tspan = (0.0, 1.0)
p = -1.0
ode = ODEProblem(f, u0, tspan, p)

param_axis = @lens _
prob = BifurcationProblem(ode, param_axis, (-1.0, 1.0))


using StaticArrays: @SVector

function closest_analytic(u, p)
    if p < 0
        return [abs(u), p]
    else
        analytic = @SVector [-sqrt(p), 0, sqrt(p)]
        _, i = findmin(abs.(u .- analytic))
        return [analytic[i], p]
    end
end

end  # module
