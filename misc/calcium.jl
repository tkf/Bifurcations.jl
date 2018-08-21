using Bifurcations
using Bifurcations.Examples: Calcium

solver = init(Calcium.prob)
@time solve!(solver)
sol = solver.sol
display(sol)

using Plots
plot = Bifurcations.plot  # workaround
plt1 = plot(sol)
display(plt1)

using Bifurcations: BifurcationProblem, special_intervals
using Setfield: @lens

point_list = sort!(special_intervals(solver), by=p->p.u0[end])
point = point_list[1]
@show point

codim2_prob = BifurcationProblem(
    point,
    solver,
    (@lens _.gca),
    (0.0, 8.0),
)
codim2_solver = init(codim2_prob)
solve!(codim2_solver)

plt2 = plot(codim2_solver.sol)
display(plt2)
