module Bautin

using DiffEqBase: ODEProblem
using Parameters: @with_kw, @unpack
using StaticArrays: SVector
using Setfield: @lens

using ...Bifurcations: BifurcationProblem

@with_kw struct BautinParam{B1, B2, S}
    β₁::B1 = -1.0
    β₂::B2 = -1.0
    s::S = -1
end

function f(u::SVector, p, t)
    @unpack β₁, β₂, s = p
    x, y = u
    #=
    ρ² = x^2 + y^2
    ρ⁴ = ρ²^2
    =#
    ρ² = x * x + y * y
    ρ⁴ = ρ² * ρ²
    a = β₂ * ρ² + s * ρ⁴
    dx = β₁ * x -      y + a * x
    dy =      x + β₁ * y + a * y
    return SVector(dx, dy)
end

# TODO: make complex Bautin normal form work
#=
function f(u::SVector, p, t)
    dz = f(u[1] + u[2] * im, p, t)
    return SVector(real(dz), imag(dz))
end

function f(z::Complex, p, t)
    @unpack β₁, β₂, s = p
    ρ² = abs2(z)
    ρ⁴ = ρ²^2
    return (β₁ + 1im) * z + β₂ * z * ρ² + s * z * ρ⁴
end
=#

f(u::A, p, t) where {A <: AbstractVector} = A(f(SVector{2}(u), p, t))

function f(du, u, p, t)
    du .= f(SVector{2}(u), p, t)
    nothing
end

make_prob(
        p = BautinParam();
        u0 = SVector(0.0, 0.0),
        tspan = (0.0, 30.0),
        ode = ODEProblem{!(u0 isa SVector)}(f, u0, tspan, p),
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
