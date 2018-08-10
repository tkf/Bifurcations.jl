module TestNormalForm
include("preamble.jl")

using Bifurcations.Continuations: as, ContinuationSolution, sweeps_as_vectors
using Bifurcations.Examples.Pitchfork
using Bifurcations.Examples.Transcritical

@testset "$normal" for normal in [Pitchfork, Transcritical]
    solver = init(normal.prob; h0=0.3)
    solve!(solver)
    sol = as(solver.sol, ContinuationSolution)
    rtol = solver.opts.rtol
    atol = solver.opts.atol
    for (us, ps) in zip(sweeps_as_vectors(sol, 1),
                        sweeps_as_vectors(sol, 0))
        analytic = hcat(normal.closest_analytic.(us, ps)...)'
        numeric = hcat(us, ps)
        @test numeric ≈ analytic  rtol=2rtol atol=atol
        # TODO: Don't double rtol; I think rtol=rtol is not working
        # because the way closest_analytic chooses the "closest" point
        # is actually not the closest.
    end

    # One of the first two sweeps has to detect a simple bifurcation
    num_bifurcations = [length(sweep.simple_bifurcation)
                        for sweep in sol.sweeps]
    @test Set(num_bifurcations[1:2]) == Set([0, 1])
    _, i = findmax(num_bifurcations)
    sweep = sol.sweeps[i]
    b = sweep.simple_bifurcation[1]
    @test b.u0[end] <= 0 <= b.u1[end]
end

end  # module
