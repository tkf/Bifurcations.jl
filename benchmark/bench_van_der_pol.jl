using BenchmarkTools

using DiffEqBase: remake
using Setfield: @lens

using Bifurcations
using Bifurcations: LimitCycleProblem
using Bifurcations.Continuations: pre_solve!, sweep!
using Bifurcations.Examples: DuffingVanDerPol

function make_van_der_pol_lc_solver()
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
        ts = ((1:n) - 1) * dt
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
        t_domain = (0.01, 1.5),
        t0 = get(DuffingVanDerPol.param_axis, ode.p),
    )

    @assert size(prob.xs0) == (2, num_mesh * degree)

    return init(
        prob;
        start_from_nearest_root = true,
        max_samples = 3,
    )
end

suite = BenchmarkGroup()
suite["sweep!(lc_solver)"] = @benchmarkable(
    sweep!(lc_solver),
    setup = (lc_solver = pre_solve!(make_van_der_pol_lc_solver())))
suite
