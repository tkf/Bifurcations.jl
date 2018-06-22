# Calcium channel model

Calcium channel model taken from PyDSTool example.  See:

* [pydstool/Tutorial_Calcium.py at master · robclewley/pydstool](https://github.com/robclewley/pydstool/blob/master/examples/Tutorial_Calcium.py)
* [Tutorial - PyDSTool Wiki](http://www2.gsu.edu/~matrhc/Tutorial.html)
* [Bifurcation Analysis · DifferentialEquations.jl](http://docs.juliadiffeq.org/latest/analysis/bifurcation.html)

Use [`QuickTypes.@qstruct_fp`][QuickTypes] to define model parameter:

```@example calcium
using QuickTypes: @qstruct_fp

@qstruct_fp CalciumParam(
    vl = -60,
    vca = 120,
    i = -220.0,
    gl = 2,
    gca = 4,
    c = 20,
    v1 = -1.2,
    v2 = 18,
)

nothing # hide
```

Define the model as in [DifferentialEquations.jl][ODEProblem]:

```@example calcium
using Parameters: @unpack

function f(u, p::CalciumParam, t)
    @unpack vl, vca, i, gl, gca, c, v1, v2 = p
    v = u[1]
    w = u[2]
    dv = (i + gl * (vl - v) - gca * 0.5 * (1 + tanh((v-v1)/v2)) * (v-vca)) / c
    dw = v-w
    return SVector(dv, dw)
end
```

Create an [`ODEProblem`][ODEProblem]:

```@example calcium
using DiffEqBase: ODEProblem
using StaticArrays: SVector

u0 = SVector(-170.0, -170.0)
tspan = (0.0, 30.0)  # ignored by Bifurcations.jl
p = CalciumParam()
ode = ODEProblem(f, u0, tspan, p)
```

Create a bifurcation problem:

```@example calcium
using Bifurcations: FixedPointBifurcationProblem
using Setfield: @lens

param_axis = @lens _.i
prob = FixedPointBifurcationProblem(ode, param_axis, (-300.0, 100.0))
nothing # hide
```

Solve it:

```@example calcium
using DiffEqBase: init, solve!

solver = init(prob)
solve!(solver)
sol = solver.sol
```

Plot it:

```@example calcium
using Plots
using Bifurcations: plot  # a workaround

plt = plot(sol)
savefig(plt, "calcium-1.png"); nothing # hide
```

![](calcium-1.png)

Find the left Saddle-Node bifurcation point:

```@example calcium
using Bifurcations: special_points

point_list = sort!(special_points(solver), by=p->p.u0[end])
point = point_list[1]
```

Numerical continuation of the Saddle-Node bifurcation point:

```@example calcium
using Bifurcations: BifurcationProblem

sn_prob = BifurcationProblem(
    point,
    solver,
    (@lens _.gca),
    (0.0, 8.0),
)
sn_solver = init(sn_prob)
solve!(sn_solver)
```

Plot the phase diagram:

```@example calcium
plt2 = plot(sn_solver.sol; vars=(5, 6))
savefig(plt2, "calcium-2.png"); nothing # hide
```

![](calcium-2.png)


[QuickTypes]: https://github.com/cstjean/QuickTypes.jl

[ODEProblem]: http://docs.juliadiffeq.org/latest/tutorials/ode_example.html
