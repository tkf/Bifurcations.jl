Pkg.init()

function required_packages(io::IO)
    filter(
        x -> x != "julia",
        [split(line)[1] for line in split(strip(readstring(io)), '\n')])
end
required_packages(path::AbstractString) = open(required_packages, path)

packages = vcat(
    required_packages("REQUIRE"),
    required_packages("test/REQUIRE"),
    ["Coverage"],
)
if get(ENV, "CI_GROUP", "") == "docs"
    append!(packages, ["Documenter", "QuickTypes"])
end
info("Installing: ", join(packages, ", "))
open(Pkg.dir("REQUIRE"), "a") do io
    writedlm(io, packages)
end

Pkg.resolve()
