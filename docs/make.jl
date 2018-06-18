using Documenter, Bifurcations

makedocs()

if get(ENV, "TRAVIS", "") == "true"
    deploydocs(
        deps   = Deps.pip("mkdocs", "python-markdown-math"),
        repo   = "github.com/tkf/Bifurcations.jl.git",
        julia  = "0.6",
    )
end
