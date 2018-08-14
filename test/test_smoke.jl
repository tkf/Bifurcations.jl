module TestSmoke
include("preamble_plots.jl")

using Bifurcations.Codim1: resolved_points, SpecialPoint
using Bifurcations: examples

@testset "smoke $name" for (name, ex) in examples()
    prob = ex.prob
    solver = init(prob)
    solve!(solver)
    points = resolved_points(solver)
    @test all(isa.(points, SpecialPoint))

    @testset "show" begin
        smoke_test_solver_show(solver)
    end

    @testset "plot" begin
        sol = solve(prob)
        @test_nothrow nullshow(plot(sol))
        @test_nothrow nullshow(plot(sol; include_points=true))
        @test_nothrow nullshow(plot(sol; bif_style=Dict()))

        @test_nothrow nullshow(plot(solver))
        @test_nothrow nullshow(plot(solver; include_points=false))

        for p in points
            @test_nothrow nullshow(plot(p))
        end
    end
end

end  # module
