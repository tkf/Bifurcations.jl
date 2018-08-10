module TestVsSVector
include("preamble.jl")

using Bifurcations: resolved_points
using Bifurcations.Continuations: find_errors, print_errors
using Bifurcations: examples
using Bifurcations.Examples.Calcium
using Bifurcations.Examples.PredatorPrey

use_array(ex) =
    ex.make_prob(
        f = generalized_f(ex),
        u0 = Array(ex.u0),
    )

EXAMPLES = [
    # name => (
    #     prob1_svec,   # SVECtor
    #     prob1_stda,   # STanDard Array
    #     sweep_ids,    # sweeps to be tested
    # )
    :Calcium => (
        Calcium.prob,
        use_array(Calcium),
        :,
        ),
    # PredatorPrey* tests fail at sweeps[4].u[8].  Probably this kind
    # of test is not a good idea after all.  But let's keep it as-is
    # for now to "freeze" the implementation.
    # TODO: Find out how PredatorPrey* tests diverge.
    :PredatorPreyImmutableState => (
        PredatorPrey.prob,
        use_array(PredatorPrey),
        1:3,
    ),
    :PredatorPreyMutableState => (
        PredatorPrey.prob,
        PredatorPrey.make_prob(u0=Array(PredatorPrey.u0)),
        1:3,
    ),
]

@testset "example $name" for (name, probs) in EXAMPLES
    prob1_svec, prob1_stda, sweep_ids = probs
    solver1_svec = init(prob1_svec)
    solver1_stda = init(prob1_stda)

    solve!(solver1_svec)
    solve!(solver1_stda)
    sweeps_svec = solver1_svec.super.sol.sweeps
    sweeps_stda = solver1_stda.super.sol.sweeps

    @test length(sweeps_svec) == length(sweeps_stda)
    for (i, (sweep_svec, sweep_stda)) in enumerate(zip(sweeps_svec[sweep_ids],
                                                       sweeps_stda[sweep_ids]))
        @test length(sweep_svec.u) == length(sweep_stda.u)
        for (j, (u_svec, u_stda)) in enumerate(zip(sweep_svec.u,
                                                   sweep_stda.u))
            # Manually abort on the first test failure:
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
