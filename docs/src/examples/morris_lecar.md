# Modified Morris-Lecar model

Modified Morris-Lecar model from [Dhooge, Govaerts, Kuznetsov (2003)]:

* [Dhooge, Govaerts, Kuznetsov (2003)].
  Numerical Continuation of Fold Bifurcations of Limit Cycles in MATCONT

[Dhooge, Govaerts, Kuznetsov (2003)]: https://doi.org/10.1007/3-540-44860-8_72

```@example morris_lecar
using Bifurcations
using Bifurcations: special_points
using Bifurcations.Codim1
using Bifurcations.Codim2
using Bifurcations.Codim2LimitCycle: FoldLimitCycleProblem
using Bifurcations.Examples: MorrisLecar

using Setfield: @lens
using Plots
```

Solve continuation of the equilibrium point:

```@example morris_lecar
solver = init(
    MorrisLecar.make_prob();
    start_from_nearest_root = true,
    max_branches = 0,
    nominal_angle_rad = 2π * (5 / 360),
)
@time solve!(solver)
```

Plot equilibriums in ``(u_1, y)``-space:

```@example morris_lecar
plt1 = plot(solver)
savefig(plt1, "morris_lecar-1.png"); nothing # hide
```

![](morris_lecar-1.png)

### Start continuation of Hopf bifurcation

```@example morris_lecar
hopf_point, = special_points(solver, Codim1.PointTypes.hopf)
```

Solve continuation of the Hopf point:

```@example morris_lecar
codim2_prob = BifurcationProblem(
    hopf_point,
    solver,
    (@lens _.z),
    (-1.0, 1.0),
)
hopf_solver1 = init(
    codim2_prob;
    nominal_angle_rad = 0.01,
)
@time solve!(hopf_solver1)
```

### Start continuation of fold bifurcation of limit cycle at Bautin bifurcation

```@example morris_lecar
bautin_point, = special_points(hopf_solver1, Codim2.PointTypes.bautin)
```

Construct a problem for fold bifurcation of the limit cycle starting
at `bautin_point`:

```@example morris_lecar
flc_prob = FoldLimitCycleProblem(
    bautin_point,
    hopf_solver1;
    period_bound = (0.0, 14.0),  # see below
    num_mesh = 120,
    degree = 4,
)
flc_solver = init(
    flc_prob;
    start_from_nearest_root = true,
    max_branches = 0,
    bidirectional_first_sweep = false,
    nominal_angle_rad = 2π * (5 / 360),
    max_samples = 500,
)
@time solve!(flc_solver)
```

Plot the limit cycles at fold bifurcation boundaries:

```@example morris_lecar
plt_state_space = plot_state_space(flc_solver)
savefig(plt_state_space, "morris_lecar-state_space.png"); nothing # hide
```

![](morris_lecar-state_space.png)

The continuation was configured to stop just before the period is
about to diverge.  Note that stopping at larger period requires larger
mesh size.

```@example morris_lecar
d_lc = length(flc_solver.prob.xs0)
plt_periods = plot(
    [u[end-1] for u in flc_solver.sol.sweeps[1].super.u],
    [u[d_lc + 1] for u in flc_solver.sol.sweeps[1].super.u],
    xlabel = "p1",
    ylabel = "period")
savefig(plt_periods, "morris_lecar-periods.png"); nothing # hide
```

![](morris_lecar-periods.png)

### Start continuation of Saddle-Node bifurcation

```@example morris_lecar
sn_point, = special_points(solver, Codim1.PointTypes.saddle_node)
```

Going back to the original continuation of the equilibrium, let's
start continuation of one of the saddle-node bifurcation:

```@example morris_lecar
sn_prob = BifurcationProblem(
    sn_point,
    solver,
    (@lens _.z),
    (-1.0, 1.0),
)
sn_solver = init(
    sn_prob;
    nominal_angle_rad = 0.01,
    max_samples = 1000,
    start_from_nearest_root = true,
)
@time solve!(sn_solver)
```

### Switching to continuation of Hopf bifurcation at Bogdanov-Takens bifurcation

```@example morris_lecar
hopf_prob2 = BifurcationProblem(
    special_points(sn_solver, Codim2.PointTypes.bogdanov_takens)[1],
    sn_solver,
)
hopf_solver2 = init(hopf_prob2)
@time solve!(hopf_solver2)
```

### Phase diagram

```@example morris_lecar
plt2 = plot()
for s in [hopf_solver1, flc_solver, sn_solver, hopf_solver2]
    n = length(s.sol.sweeps[1].super.u[1])
    plot!(plt2, s, vars=(n - 1, n))
end
plot!(plt2, ylim=(0.03, 0.15), xlim=(-0.05, 0.2))

savefig(plt2, "morris_lecar-2.png"); nothing # hide
```

![](morris_lecar-2.png)
