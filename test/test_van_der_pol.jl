module TestVanDerPol
include("preamble.jl")

using Bifurcations: LimitCycleProblem
using Bifurcations.Examples: VanDerPol
using Bifurcations.Continuations: as, ContinuationSolution, sweeps_as_vectors

# Create a limit cycle solution.
using OrdinaryDiffEq: Tsit5
using DiffEqBase: remake
param = VanDerPol.VanDerPolParam(
    d = 0.1,
)
ode = remake(
    VanDerPol.ode,
    u0 = [1.0, 1.8],
    p = param,
    tspan = (0.0, 100),
)
sol = solve(ode, Tsit5())

using Roots: find_zero
# I know that t1 - t0 ~ 6.2 so I can do this:
t0 = find_zero((t) -> sol(t)[1] - 1, (80, 83))
t1 = find_zero((t) -> sol(t)[1] - 1, (t0 + 3, t0 + 7))

x0 = sol(t0)
@assert all(isapprox.(x0, sol(t1); rtol=1e-2))

num_mesh = 20
degree = 3
prob = LimitCycleProblem(
    ode, VanDerPol.param_axis, VanDerPol.t_domain,
    num_mesh, degree;
    x0 = x0,
    l0 = t1 - t0,
    t_domain = (0.01, 1.5), # bound it so that it works with above `num_mesh`
    de_args = [Tsit5()],
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
