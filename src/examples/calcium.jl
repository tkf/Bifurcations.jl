module Calcium

using DiffEqBase: ODEProblem
using Parameters: @with_kw, @unpack
using StaticArrays: SVector
using Setfield: @lens
import Setfield

using ...Bifurcations: FixedPointBifurcationProblem

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
function f(u, p::CalciumParam, t)
    @unpack vl, vca, i, gl, gca, c, v1, v2 = p
    v = u[1]
    w = u[2]
    dv = (i + gl * (vl - v) - gca * 0.5 * (1 + tanh((v-v1)/v2)) * (v-vca)) / c
    dw = v-w
    return SVector(dv, dw)
end

u0 = SVector(-170.0, -170.0)
tspan = (0.0, 30.0)
p = CalciumParam()
ode = ODEProblem(f, u0, tspan, p)

param_axis = @lens _.i
prob = FixedPointBifurcationProblem(ode, param_axis, (-300.0, 100.0))

end  # module
