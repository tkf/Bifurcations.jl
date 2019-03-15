using Documenter

# https://docs.travis-ci.com/user/environment-variables/#Default-Environment-Variables
if get(ENV, "TRAVIS", "") != "true"
    # Don't do anything outside Travis CI
elseif get(ENV, "CI_GROUP", "") != "docs"
    @info("Skipping deploy since CI_GROUP != docs.")
elseif startswith(get(ENV, "TRAVIS_BRANCH", ""), "pre/")
    # For branches pre/*, deploy them into gh-pages.pre.
    branch = ENV["TRAVIS_BRANCH"]
    deploydocs(
        deps   = Deps.pip("mkdocs", "python-markdown-math"),
        make   = () -> run(setenv(`mkdocs build`, dir=@__DIR__)),
        repo   = "github.com/tkf/Bifurcations.jl.git",
        branch = "gh-pages.pre",
        devbranch = branch,
        target = "site",
        root   = @__DIR__,
    )
else
    deploydocs(
        deps   = Deps.pip("mkdocs", "python-markdown-math"),
        make   = () -> run(setenv(`mkdocs build`, dir=@__DIR__)),
        repo   = "github.com/tkf/Bifurcations.jl.git",
        target = "site",
        root   = @__DIR__,
    )
end
