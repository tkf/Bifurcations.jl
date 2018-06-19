module TestSmoke
using Compat
using Base.Test
using Plots
using Bifurcations
using Bifurcations.Examples: PROBLEMS
include("utils.jl")

function nullshow(plt)
    show(devnull, MIME("image/png"), plt)
end

for prob in PROBLEMS
    @test_nothrow begin
        sol = solve(prob)
        plt = plot(sol)
        nullshow(plt)
    end
end

end  # module
