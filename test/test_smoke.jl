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

@testset "solving!" begin
    (_, ex), = examples()
    prob = ex.prob
    solver = init(prob)
    points = []
    solving!(solver) do x
        push!(points, (; map(n -> n => getproperty(x, n), propertynames(x))...))
    end
    @test length(points) > 0
    @test keys(points[1]) == (:i_sweep, :i_point, :u, :direction, :simple_bifurcation)
    @test all(p.u isa AbstractVector for p in points)
end

@testset "include(gpu.jl)" begin
    m = Module()
    @test begin
        Base.include(m, "../examples/gpu.jl")
        m.solver isa Bifurcations.Continuations.ContinuationSolver
    end
end

end  # module
