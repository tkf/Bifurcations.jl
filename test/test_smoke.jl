module TestSmoke
using Base.Test
using Plots
using Bifurcations
using Bifurcations.Examples: PROBLEMS
include("utils.jl")

for prob in PROBLEMS
    @test_nothrow begin
        sol = solve(prob)
        plot(sol)
    end
end

end  # module
