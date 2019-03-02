# Run some non-test portion of standard script here.  See also:
# https://github.com/travis-ci/travis-build/blob/master/lib/travis/build/script/julia.rb

using Pkg

# Manually installing other CI dependencies
packages = []
if get(ENV, "TRAVIS", "") == "true"
    append!(packages, ["Coverage"])
end
if get(ENV, "CI_GROUP", "") == "docs"
    append!(packages, [
        "Documenter", "DocumenterMarkdown", "QuickTypes", "Roots",
    ])
end
specs = [
    PackageSpec.(packages)...
]
@info string("Installing:\n", join(specs, "\n"))
Pkg.add(specs)
