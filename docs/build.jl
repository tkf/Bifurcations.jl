printstyled(let
    message = " Warnings below (if any) are fine. "
    margin = (displaysize(stdout)[2] - length(message)) ÷ 2
    ("=" ^ margin) * message * ("=" ^ margin)
end, color=:blue)
println()
flush(stdout)
using Documenter
import DocumenterMarkdown
using DiffEqBase: AbstractODEProblem
using Setfield: Lens
import Plots
import OrdinaryDiffEq
flush(stdout)
flush(stderr)
printstyled("=" ^ displaysize(stdout)[2], color=:blue)
println()

using Bifurcations

let strict = get(ENV, "DOCS_STRICT", "yes") == "yes"
    @info """
    `makedocs` with:
    strict = $strict
    """
    makedocs(
        strict = strict,
        format = DocumenterMarkdown.Markdown()
    )
end
