module TestSmoke
using Base.Test

import GR  # Workaround Plots.jl world age problem
using Plots

using Bifurcations
using Bifurcations.Codim1: resolved_points, SpecialPoint
using Bifurcations.Examples: PROBLEMS
include("utils.jl")

function nullshow(plt::Plots.Plot)
    nullshow(MIME("image/png"), plt)
end

@testset "smoke PROBLEMS[$i]" for (i, prob) in enumerate(PROBLEMS)
    solver = init(prob)
    solve!(solver)
    points = resolved_points(solver)
    @test all(isa.(points, SpecialPoint))

    @testset "show" begin
        @test_nothrow nullshow(solver)
        @test_nothrow nullshow(solver.sol)

        for sweep in solver.sol.sweeps
            @test_nothrow nullshow(sweep)
            for point in sweep.special_points
                @test_nothrow nullshow(point)
            end
        end
    end

    @testset "plot" begin
        sol = solve(prob)
        @test_nothrow nullshow(plot(sol))
        @test_nothrow nullshow(plot(sol; include_points=true))

        @test_nothrow nullshow(plot(solver))
        @test_nothrow nullshow(plot(solver; include_points=true))

        for p in points
            @test_nothrow nullshow(plot(p))
        end
    end
end

end  # module
