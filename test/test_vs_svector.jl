module TestVsSVector
include("preamble.jl")

using DiffEqBase: remake

using Bifurcations: Codim1, resolved_points
using Bifurcations.Continuations: find_errors, print_errors
using Bifurcations: examples
using Bifurcations.Examples: Calcium

EXAMPLES = [
    :Calcium => Calcium,
]

@testset "example $name" for (name, ex) in EXAMPLES
    prob1_svec = ex.prob
    prob1_stda = ex.make_prob(  # STanDard Array
        ode = remake(ex.ode;
                     f = generalized_f(ex),
                     u0 = Array(ex.u0)),
    )

    solver1_svec = init(prob1_svec)
    solver1_stda = init(prob1_stda)

    solve!(solver1_svec)
    solve!(solver1_stda)
    sweeps_svec = solver1_svec.super.sol.sweeps
    sweeps_stda = solver1_stda.super.sol.sweeps

    @test length(sweeps_svec) == length(sweeps_stda)
    for (i, (sweep_svec, sweep_stda)) in enumerate(zip(sweeps_svec,
                                                       sweeps_stda))
        @test length(sweep_svec.u) == length(sweep_stda.u)
        for (j, (u_svec, u_stda)) in enumerate(zip(sweep_svec.u,
                                                   sweep_stda.u))
            if ! (u_svec ≈ u_stda)
                println("sweeps[$i].u[$j]")
                @test u_svec ≈ u_stda
                break
            else
                @test u_svec ≈ u_stda
            end
        end
    end
end

end  # module
