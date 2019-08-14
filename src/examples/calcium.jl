module Calcium

using DiffEqBase: ODEProblem, ODEFunction
using ForwardDiff
using Parameters: @with_kw, @unpack
using StaticArrays: SVector
using Setfield: @lens
import Setfield

using ...Bifurcations: BifurcationProblem
using ...Codim2: DiffEqCodim2Problem

@with_kw struct CalciumParam{
        vlType, vcaType, iType, glType, gcaType, cType, v1Type, v2Type}
    vl::vlType = -60
    vca::vcaType = 120
    i::iType = -220.0
    gl::glType = 2
    gca::gcaType = 4
    c::cType = 20
    v1::v1Type = -1.2
    v2::v2Type = 18
end

# http://www2.gsu.edu/~matrhc/Tutorial.html
# https://github.com/robclewley/pydstool/blob/master/examples/Tutorial_Calcium.py
# http://docs.juliadiffeq.org/latest/analysis/bifurcation.html
function f(u::SVector, p::CalciumParam, t)
    @unpack vl, vca, i, gl, gca, c, v1, v2 = p
    v = u[1]
    w = u[2]
    dv = (i + gl * (vl - v) - gca * 0.5 * (1 + tanh((v-v1)/v2)) * (v-vca)) / c
    dw = v-w
    return SVector(dv, dw)
end

function f(du, u, p, t)
    du .= f(SVector{2}(u), p, t)
    nothing
end

function jac(J, u, p, t)
    J .= ForwardDiff.jacobian(x -> f(x, p, t), SVector{2}(u))
    nothing
end

const _f = f
const _jac = jac


make_prob(
        p = CalciumParam();
        f = _f,
        u0 = SVector(-170.0, -170.0),
        tspan = (0.0, 30.0),
        ode = ODEProblem{!(u0 isa SVector)}(f, u0, tspan, p),
        param_axis = (@lens _.i),
        t_domain = (-300.0, 100.0),
        kwargs...) =
    BifurcationProblem(ode, param_axis, t_domain;
                       kwargs...)

prob = make_prob()
ode = prob.p.de_prob
u0 = ode.u0
param_axis = prob.p.param_axis


function make_codim2_prob(
        p = CalciumParam();
        u0 = SVector(-170.0, -170.0),
        tspan = (0.0, 30.0),
        ode = ODEProblem{!(u0 isa SVector)}(f, u0, tspan, p),
        param_axis1 = (@lens _.i),
        param_axis2 = (@lens _.gca),
        t_domain = ([-300.0, 0.0], [100.0, 8.0]),
        kwargs...)
    return DiffEqCodim2Problem(
        ode,
        param_axis1,
        param_axis2,
        t_domain;
        kwargs...)
end

make_ode_jac(
    p = CalciumParam();
    u0 = [-170.0, -170.0],
    tspan = (0.0, 30.0),
    f = _f,
    jac = _jac,
) =
    ODEProblem(ODEFunction{true}(f; jac = jac), u0, tspan, p)

end  # module
