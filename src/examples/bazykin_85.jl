module Bazykin85

using DiffEqBase: ODEProblem
using Parameters: @with_kw, @unpack
using StaticArrays: SVector
using Setfield: @lens
import Setfield

using ...Bifurcations: BifurcationProblem

@with_kw struct Bazykin85Param{A, E, G, D}
    α::A = 0.1
    ϵ::E = 0.01  # ≪ 1 ?
    γ::G = 1.0
    δ::D = 0.1
end


function f(x, p, t)
    @unpack α, ϵ, γ, δ = p
    return SVector(
             x[1] - x[1] * x[2] / (1 + α * x[1]) - ϵ * x[1]^2,
        -γ * x[2] + x[1] * x[2] / (1 + α * x[1]) - δ * x[2]^2,
    )
end

tspan = (0.0, 30.0)
p = Bazykin85Param()
# u0 = SVector(0.0, 0.0)
u0 = SVector(1 / p.ϵ, 0.0)
ode = ODEProblem(f, u0, tspan, p)

param_axis = @lens _.α
prob = BifurcationProblem(ode, param_axis, (0.01, 1.5),
                          phase_space = (SVector(0.0, 0.0),  # u_min
                                         SVector(Inf, Inf)))   # u_max

end  # module
