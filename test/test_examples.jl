module TestExamples
include("preamble.jl")

using Bifurcations: Codim1, resolved_points
using Bifurcations.Continuations: find_errors, print_errors
using Bifurcations: examples
using Bifurcations.Examples.Calcium: CalciumParam
using Bifurcations.Examples.Bazykin85: Bazykin85Param

make_codim2(::Real, point, solver1) = nothing
make_codim2(::Tuple, point, solver1) = nothing  # PredatorPrey

make_codim2(::CalciumParam, point, solver1) =
    BifurcationProblem(
        point,
        solver1,
        (@lens _.gca),
        (0.0, 8.0),
    )

function make_codim2(::Bazykin85Param, point, solver1)
    if (point.point_type == Codim1.PointTypes.saddle_node ||
        isapprox(point.u[2], 0; atol=1e-3))
        return nothing
    end
    return BifurcationProblem(
        point,
        solver1,
        (@lens _.δ),
        (0.0, 10.0),
    )
end

@testset "example $name" for (name, ex) in examples()
    prob1 = ex.prob
    solver1 = init(prob1)
    solve!(solver1)
    errors = find_errors(solver1)
    print_errors(errors)
    @test isempty(errors)

    for point in resolved_points(solver1)
        prob2 = make_codim2(prob1.p.de_prob.p, point, solver1)
        if prob2 === nothing
            continue
        end
        solver2 = init(prob2)
        solve!(solver2)
        errors = find_errors(solver2)
        print_errors(errors)
        @test isempty(errors)
    end
end

end  # module
