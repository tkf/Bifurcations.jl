module TestExamples
include("preamble.jl")

using DiffEqBase: remake

using Bifurcations: Codim1, resolved_points
using Bifurcations.Codim2: BackReferencingAS, NormalizingAS
using Bifurcations.Continuations: find_errors, print_errors
using Bifurcations: examples
using Bifurcations.Examples: Calcium
using Bifurcations.Examples.Calcium: CalciumParam
using Bifurcations.Examples.Bazykin85: Bazykin85Param
using Bifurcations.Examples.DuffingVanDerPol: DuffingVanDerPolParam
using Bifurcations.Examples.Bautin: BautinParam
using Bifurcations.Examples.Cusp: CuspParam

struct Example{P}
    prob::P
end

make_codim2(::Real, point, solver1) = []
make_codim2(::Tuple, point, solver1) = []  # PredatorPrey
make_codim2(::DuffingVanDerPolParam, point, solver1) = []
make_codim2(::CuspParam, point, solver1) = []  # TODO: define something

make_codim2(::CalciumParam, point, solver1) = [
    BifurcationProblem(
        point,
        solver1,
        (@lens _.gca),
        (0.0, 8.0),
    )
]

function make_codim2(::Bazykin85Param, point, solver1)
    if (point.point_type == Codim1.PointTypes.saddle_node ||
        isapprox(point.u[2], 0; atol=1e-3))
        return []
    end
    return [BifurcationProblem(
        point,
        solver1,
        (@lens _.δ),
        (0.0, 10.0),
    )]
end

function make_codim2(::BautinParam, point, solver1)
    args = (
        point,
        solver1,
        (@lens _.β₂),
        (-2.0, 2.0),
    )
    return [
        BifurcationProblem(
            args...;
            augmented_system = BackReferencingAS(),
        ),
        BifurcationProblem(
            args...;
            augmented_system = NormalizingAS(),
        ),
    ]
end

EXAMPLES = [
    examples()...,
    :CalciumSTDArray => Example(Calcium.make_prob(
        ode = remake(Calcium.ode;
                     f = generalized_f(Calcium),
                     u0 = Array(Calcium.u0)),
    )),
]

@testset "example $name" for (name, ex) in EXAMPLES
    prob1 = ex.prob
    solver1 = init(prob1)
    solve!(solver1)
    errors = find_errors(solver1)
    print_errors(errors)
    @test isempty(errors)

    @testset "$i-th $(point.point_type); $j-th problem" for
            (i, point) in enumerate(resolved_points(solver1)),
            (j, prob2) in enumerate(
                make_codim2(prob1.p.de_prob.p, point, solver1))
        solver2 = init(prob2)
        solve!(solver2)
        errors = find_errors(solver2)
        print_errors(errors)
        @test isempty(errors)
    end
end

end  # module
