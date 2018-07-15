using ..Continuations: ContinuationOptions
using ..Codim1
using ..Codim1: Codim1Solver, resolve_point
using ..BifurcationsBase: AbstractSpecialPoint, SpecialPoint
using ...Bifurcations: FixedPointBifurcationProblem, DiffEqWrapper
using DiffEqBase: ODEProblem, remake, solve

"""
    LimitCycleProblem(point::AbstractSpecialPoint,
                      solver)

Construct limit cycle continuation (codimension-1 bifurcation) problem
given a bifurcation `point`.
"""
function LimitCycleProblem(point::AbstractSpecialPoint,
                           solver::Codim1Solver;
                           kwargs...)

    cd1_prob = solver.prob :: FixedPointBifurcationProblem
    de_wrapper = cd1_prob.p :: DiffEqWrapper
    de_prob = de_wrapper.de_prob :: ODEProblem

    return LimitCycleProblem(
        resolve_point(point, solver),
        solver.opts;
        de_prob = de_prob,
        param_axis = de_wrapper.param_axis,
        kwargs...)
end

function LimitCycleProblem(point::SpecialPoint{Continuous},
                           opts::ContinuationOptions = ContinuationOptions();
                           num_mesh::Int = error("num_mesh required"),
                           degree::Int = error("degree required"),
                           amp = opts.h0,
                           kwargs...)
    @assert point.point_type == Codim1.PointTypes.hopf

    center = point.u[1:end - 1]
    t0 = point.u[end]

    vals, vecs = eig(Codim1.ds_jacobian(point.J))
    _, idx0 = findmin(abs.(real.(vals)))  # real closest to zero
    v0 = normalize(vecs[:, idx0])
    w0 = imag(vals[idx0])
    if w0 < 0
        w0 = -w0
        v0 = conj(v0)
    end
    @assert w0 > opts.atol

    # TODO: use predictor-corrector rather than guessing amplitude
    vs = real(v0) .* amp
    vc = imag(v0) .* amp

    points = num_mesh * degree
    rads = range(0, 2π / points, points)
    xs0 = center .+ (vs * sin.(rads)' .+
                     vc * cos.(rads)')

    l0 = 2π / w0  # period
    return LimitCycleProblem(;
        xs0 = xs0,
        l0 = l0,
        t0 = t0,
        num_mesh = num_mesh,
        degree = degree,
        kwargs...)
end

"""
    LimitCycleProblem(ode, param_axis, t_domain,
                      num_mesh, degree;
                      <keyword arguments>)

Solve `ode` and use it's solution as the initial limit cycle of the
continuation problem.

# Arguments
- `ode::ODEProblem`
- `param_axis::Lens`
- `t_domain::Tuple`
- `num_mesh::Integer`
- `degree::Integer`
"""
function LimitCycleProblem(de_prob::ODEProblem,
                           param_axis::Lens,
                           t_domain::Tuple,
                           num_mesh::Int,
                           degree::Int;
                           x0 = de_prob.u0,
                           l0 = de_prob.tspan[end] - de_prob.tspan[1],
                           time_offset = de_prob.tspan[1],
                           de_args = [],
                           de_opts = [],
                           kwargs...)
    ode = remake(
        de_prob;
        u0 = x0,
        tspan = (time_offset, time_offset + l0),
    )
    sol = solve(ode, de_args...; de_opts...)

    xs0 = similar(ode.u0, (length(ode.u0), num_mesh * degree))
    for (i, t) in enumerate(linspace(0, l0, num_mesh * degree))
        xs0[:, i] = sol(t)
    end

    return LimitCycleProblem(;
        de_prob = ode,
        xs0 = xs0,
        l0 = l0,
        num_mesh = num_mesh,
        degree = degree,
        param_axis = param_axis,
        t_domain = t_domain,
        t0 = get(param_axis, ode.p),
        kwargs...)
end
