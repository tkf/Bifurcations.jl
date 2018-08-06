using Documenter

# https://docs.travis-ci.com/user/environment-variables/#Default-Environment-Variables
if get(ENV, "TRAVIS", "") != "true"
    # Don't do anything outside Travis CI
elseif get(ENV, "CI_GROUP", "") != "docs"
    info("Skipping deploy since CI_GROUP != docs.")
elseif startswith(get(ENV, "TRAVIS_BRANCH", ""), "pre/")
    # For branches pre/*, deploy them into gh-pages.pre.
    branch = ENV["TRAVIS_BRANCH"]
    deploydocs(
        deps   = Deps.pip("mkdocs", "python-markdown-math"),
        repo   = "github.com/tkf/Bifurcations.jl.git",
        julia  = "0.6",
        branch = "gh-pages.pre",
        latest = branch,
    )
else
    deploydocs(
        deps   = Deps.pip("mkdocs", "python-markdown-math"),
        repo   = "github.com/tkf/Bifurcations.jl.git",
        julia  = "0.6",
    )
end
