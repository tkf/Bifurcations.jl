print_with_color(:blue, let
    message = " Warnings below (if any) are fine. "
    margin = (displaysize(STDOUT)[2] - length(message)) ÷ 2
    ("=" ^ margin) * message * ("=" ^ margin)
end)
println()
flush(STDOUT)
using Documenter
using DiffEqBase: AbstractODEProblem
using Setfield: Lens
import Plots
import Jacobi
import OrdinaryDiffEq
flush(STDOUT)
flush(STDERR)
print_with_color(:blue, "=" ^ displaysize(STDOUT)[2])
println()

using Bifurcations
using Bifurcations: plot, plot!  # a workaround

makedocs(
    strict = true,
)
