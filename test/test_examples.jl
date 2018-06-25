module TestExamples
include("preamble.jl")

using Bifurcations.Continuations: find_errors, print_errors
using Bifurcations.Examples: PROBLEMS
using Bifurcations.Examples.Calcium: CalciumParam

make_codim2(::Real, point, solver1) = nothing
make_codim2(::Tuple, point, solver1) = nothing  # PredatorPrey

make_codim2(::CalciumParam, point, solver1) =
    BifurcationProblem(
        point,
        solver1,
        (@lens _.gca),
        (0.0, 8.0),
    )

@testset "PROBLEMS[$i]" for (i, prob1) in enumerate(PROBLEMS)
    solver1 = init(prob1)
    solve!(solver1)
    errors = find_errors(solver1)
    print_errors(errors)
    @test isempty(errors)

    for point in special_points(solver1)
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
