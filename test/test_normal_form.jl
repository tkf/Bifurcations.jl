module TestNormalForm

using Base.Test
using Bifurcations
using Bifurcations.Continuations: sweeps_as_vectors, isalmostzero
using Bifurcations.Examples: Pitchfork, Transcritical

@testset "$normal" for normal in [Pitchfork, Transcritical]
    solver = init(normal.prob; h0=0.3)
    solve!(solver)
    sol = solver.sol
    rtol = solver.opts.rtol
    atol = solver.opts.atol
    for (us, ps) in zip(sweeps_as_vectors(sol, 1),
                        sweeps_as_vectors(sol, 0))
        @test isalmostzero(normal.deviation.(us, ps), rtol, atol)
    end

    # One of the first two sweeps has to detect a simple bifurcation
    num_bifurcations = [length(sweep.simple_bifurcation)
                        for sweep in sol.sweeps]
    @test Set(num_bifurcations[1:2]) == Set([0, 1])
    _, i = findmax(num_bifurcations)
    sweep = sol.sweeps[i]
    p = [u[end] for u in sweep.u]
    b = sweep.simple_bifurcation[1]
    @test p[b] <= 0 <= p[b + 1]
end

end  # module
