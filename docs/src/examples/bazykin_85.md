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
using Bifurcations: special_points
using Setfield: @lens

point_list = sort!(special_points(solver), by=p->p.u0[end])

codim2_solvers = []
for point in point_list[2:3]
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
```

Merge two continuations and draw the bifurcation diagram:

```@example bazykin85
plt2 = plot()
for s in codim2_solvers
    n = length(s.sol.sweeps[1].super.u[1])
    plot!(plt2, s, vars=(n - 1, n))
end
savefig(plt2, "bazykin85-2.png"); nothing # hide
```

![](bazykin85-2.png)
