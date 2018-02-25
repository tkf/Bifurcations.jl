module TestSmoke
using Base.Test
using Bifurcations
using Bifurcations.Examples: PROBLEMS
include("utils.jl")

for prob in PROBLEMS
    @test_nothrow solve(prob)
end

end  # module
