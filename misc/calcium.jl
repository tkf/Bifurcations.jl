using Bifurcations
using Bifurcations.Examples: Calcium
solver = init(Calcium.prob)
@time solve!(solver)
sol = solver.sol
display(sol)
using Plots
plot(sol)
