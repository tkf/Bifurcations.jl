printstyled(let
    message = " Warnings below (if any) are fine. "
    margin = (displaysize(stdout)[2] - length(message)) ÷ 2
    ("=" ^ margin) * message * ("=" ^ margin)
end, color=:blue)
println()
flush(stdout)
using Documenter
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
        # modules = [Bifurcations],
        pages = [
            "Home" => "index.md",
            "API" => "api.md",
            "Examples" => [
                hide("examples.md"),
                "examples/calcium.md",
                "examples/van_der_pol.md",
                "examples/morris_lecar.md",
            ],
            hide("Internals" => "internals.md"),
        ],
        repo = "https://github.com/tkf/Bifurcations.jl/blob/{commit}{path}#L{line}",
        sitename = "Bifurcations.jl",
        authors = "Takafumi Arakaki",
        assets = [],
        strict = strict,
    )
end
