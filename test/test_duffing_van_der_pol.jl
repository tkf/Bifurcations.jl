module TestVanDerPol
include("preamble.jl")

using Bifurcations: LimitCycleProblem
using Bifurcations.Examples.DuffingVanDerPol
using Bifurcations.Continuations: as, ContinuationSolution, sweeps_as_vectors

using DiffEqBase: remake
param = DuffingVanDerPol.DuffingVanDerPolParam(
    d = 0.1,
)
ode = remake(
    DuffingVanDerPol.ode,
    p = param,
)

num_mesh = 20
degree = 3

# Start continuation of the limit cycle at the small damping limit.
xs0 = let
    n = num_mesh * degree
    dt = 2π / n
    ts = ((1:n) .- 1) .* dt
    hcat(
        2 .* cos.(ts),
        -2 .* sin.(ts),
    )'
end

prob = LimitCycleProblem(
    de_prob = ode;
    param_axis = DuffingVanDerPol.param_axis,
    num_mesh = num_mesh,
    degree = degree,
    xs0 = xs0,
    l0 = 2π,
    t_domain = (0.01, 1.5), # bound it so that it works with above `num_mesh`
    t0 = get(ode.p, DuffingVanDerPol.param_axis),
)

@test size(prob.xs0) == (2, num_mesh * degree)

solver = init(
    prob;
    start_from_nearest_root = true,
)
solve!(solver)

parameters, = let sol = as(solver.sol, ContinuationSolution)
    sweeps_as_vectors(sol, length(sol.sweeps[1].u[1]))
end
@test minimum(parameters) <= prob.t_domain[1]
@test maximum(parameters) >= prob.t_domain[end]

end  # module
