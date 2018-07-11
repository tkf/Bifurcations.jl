# Bazykin's predator-prey system

```@example bazykin85
using Bifurcations
using Bifurcations.Examples: Bazykin85

solver = init(Bazykin85.prob)
solve!(solver)
```

Plot codimension-one bifurcations:

```@example bazykin85
using Bifurcations: plot, plot!  # workaround
using Plots

plt1 = plot(solver.sol)
savefig(plt1, "bazykin85-1.png"); nothing # hide
```

![](bazykin85-1.png)

Let's follow the Hopf and Saddle-Node bifurcations:

```@example bazykin85
using Bifurcations: Codim1, special_points
using Setfield: @lens

sn_point, = sort!(
    special_points(solver, Codim1.PointTypes.saddle_node),
    by=p->p.u0[1])
hopf_point, = special_points(solver, Codim1.PointTypes.hopf)

point_list = [sn_point, hopf_point]

codim2_solvers = []
for point in point_list
    @show point

    codim2_prob = BifurcationProblem(
        point,
        solver,
        (@lens _.δ),
        (0.0, 10.0),
    )
    codim2_solver = init(
        codim2_prob;
        nominal_angle_rad = 0.01,
        max_samples = 1000,
    )
    push!(codim2_solvers, codim2_solver)
    solve!(codim2_solver)

    @show codim2_solver
end

sn_solver1, hopf_solver1 = codim2_solvers

nothing # hide
```

Switch to Hopf bifurcation via Bogdanov-Takens bifurcation:

```@example bazykin85
using Bifurcations: Codim2
bt_point, = sort(
    special_points(sn_solver1,
                   Codim2.PointTypes.bogdanov_takens);
    by = p -> p.u0[end - 1],  # α
    rev = true)
hopf_prob = BifurcationProblem(bt_point, sn_solver1)
hopf_solver2 = init(hopf_prob)
solve!(hopf_solver2)
@show hopf_solver2

push!(codim2_solvers, hopf_solver2)

nothing # hide
```

Merge continuations and draw the bifurcation diagram:

```@example bazykin85
plt2 = plot()
for s in codim2_solvers
    n = length(s.sol.sweeps[1].super.u[1])
    plot!(plt2, s, vars=(n - 1, n))
end
savefig(plt2, "bazykin85-2.png"); nothing # hide
```

![](bazykin85-2.png)
