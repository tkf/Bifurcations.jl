# Run some non-test portion of standard script here.  See also:
# https://github.com/travis-ci/travis-build/blob/master/lib/travis/build/script/julia.rb

using Pkg
Pkg.build()

# Manually installing other CI dependencies
packages = ["Coverage"]
if get(ENV, "CI_GROUP", "") == "docs"
    append!(packages, ["Documenter", "QuickTypes", "Roots"])
end
@info string("Installing: ", join(packages, ", "))
Pkg.add(packages)
