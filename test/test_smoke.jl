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

@testset "warn" begin
    plot = Plots.plot
    solver = init(Bifurcations.Examples.Calcium.prob)
    solve!(solver)
    sol = solver.sol

    shows = [
        () -> nullshow(plot(sol.sweeps[1]; include_points=true)),
        () -> nullshow(plot(sol; include_points=true)),
        () -> nullshow(plot(solver)),
    ]

    for show_plot in shows
        @test_logs(
            (:warn, r"include_points = true"),
            match_mode=:any,
            show_plot())
    end
end

end  # module
