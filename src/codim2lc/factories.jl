using DiffEqBase: AbstractODEProblem, remake
using StaticArrays: SVector
using Setfield: get, set

using ..Continuations: init, solve!, residual_jacobian
using ..BifurcationsBase: AbstractSpecialPoint, SpecialPoint, special_points
using ..Codim1
using ..Codim1: resolve_point
using ..Codim2
using ..Codim2: Codim2Solver, DiffEqCodim2Problem, first_lyapunov_coefficient
using ..Codim1LimitCycle: Codim1LCSolver

# See also:
# * [[../codim2/problem.jl::function BifurcationProblem]]
# * [[../codim2/switch.jl::function BifurcationProblem]]
# * [[../codim1lc/factories.jl::function LimitCycleProblem]]
function FoldLimitCycleProblem(point::AbstractSpecialPoint,
                               solver::Codim2Solver;
                               kwargs...)
    @assert point.point_type == Codim2.PointTypes.bautin
    @assert timekind(point) isa Continuous

    solver.prob :: DiffEqCodim2Problem
    de_prob = solver.prob.de_prob :: AbstractODEProblem

    # TODO: Stop solving sub-problems in the constructor:
    lc_prob = bautin_to_subcritical_lc(point, solver; kwargs...)
    lc_solver = init(
        lc_prob;
        start_from_nearest_root = true,
        max_branches = 0,
        bidirectional_first_sweep = false,
        nominal_angle_rad = 2π * (1 / 360),  # TODO: should it be the default?
    )
    solve!(lc_solver)

    flc_points = special_points(lc_solver,
                                Codim1LimitCycle.PointTypes.saddle_node)
    flc_point = resolve_point(flc_points[1], lc_solver)
    if length(flc_points) > 1
        param2 = get(solver.prob.param_axis2, lc_solver.prob.de_prob.p)
        @warn """
        $(length(flc_points)) fold bifurcations are found.
        Using the first one found:
            $(solver.prob.param_axis1) : $(flc_point.u[end])
            $(solver.prob.param_axis2) : $(param2) (fixed)
        """
    end

    n = length(lc_prob.xs0)
    vals, vecs = eig(@view flc_point.J[1:n, 1:n])
    _, idx0 = findmin(abs.(real.(vals)))  # real closest to zero
    @assert all(x -> imag(x) ≈ 0, vecs[:, idx0])
    v0 = normalize(real(@view vecs[:, idx0]))

    xs0 = reshape(flc_point.u[1:length(lc_prob.xs0)], size(lc_prob.xs0))
    l0 = flc_point.u[length(lc_prob.xs0) + 1]
    vs0 = reshape(v0, size(lc_prob.xs0))
    t0 = SVector(flc_point.u[end],
                 get(solver.prob.param_axis2, lc_solver.prob.de_prob.p))

    return FoldLimitCycleProblem(
        lc_prob;
        xs0 = xs0,
        l0 = l0,
        vs0 = vs0,
        dl0 = 0.0,  # ???
        t0 = t0,
        param_axis1 = solver.prob.param_axis1,
        param_axis2 = solver.prob.param_axis2,
        t_domain = solver.prob.t_domain,
    )
end

function BifurcationProblem(point::AbstractSpecialPoint,
                            solver::Codim1LCSolver,
                            param_axis2::Lens,
                            t2_domain::Tuple;
                            kwargs...)

    cd1_prob = solver.prob :: LimitCycleProblem
    de_prob = cd1_prob.de_prob :: AbstractODEProblem

    xs0 = reshape(point.u[1:end-2], size(cd1_prob.xs0))

    # Initialize perturbation direction
    vs0 = xs0 .- mean(xs0, 2)
    normalize!(view(vs0, :))

    l0 = point.u[end - 1]
    t0 = SVector(point.u[end], get(param_axis2, de_prob.p))
    t_domain = (
        SVector(cd1_prob.t_domain[1], t2_domain[1]),
        SVector(cd1_prob.t_domain[2], t2_domain[2]),
    )
    param_axis1 = cd1_prob.param_axis

    return FoldLimitCycleProblem(
        cd1_prob;
        xs0 = xs0,
        l0 = l0,
        vs0 = vs0,
        dl0 = 0.0,  # ???
        t0 = t0,
        param_axis1 = param_axis1,
        param_axis2 = param_axis2,
        t_domain = t_domain,
        kwargs...)
end


function bautin_to_subcritical_lc(
        point::AbstractSpecialPoint,
        solver::Codim2Solver;
        kwargs...)

    point :: SpecialPointInterval  # TODO: support SpecialPoint
    @assert point.point_type == Codim2.PointTypes.bautin
    @assert timekind(point) isa Continuous

    solver.prob :: DiffEqCodim2Problem
    de_prob = solver.prob.de_prob :: AbstractODEProblem
    prob_cache = as(solver.cache, ContinuationCache).prob_cache

    coeff_1 = first_lyapunov_coefficient(prob_cache, point.u1, point.J1)
    if coeff_1 > 0  # u1 is subcritical
        u_sc = point.u1
        J_sc = point.J1
    else
        u_sc = point.u0
        _H, J_sc = residual_jacobian(point.u0, prob_cache)
    end
    coeff_0 = first_lyapunov_coefficient(
        prob_cache, point.u0, residual_jacobian(point.u0, prob_cache)[2])
    @assert coeff_0 * coeff_1 < 0

    # [[../codim2/switch.jl::Manually cast]]
    xtype = typeof(de_prob.u0)

    N = length(de_prob.u0)
    x0 = xtype(u_sc[1:N])
    t0 = u_sc[end-1]

    hopf_point = SpecialPoint(
        timekind(point),
        Codim1.PointTypes.hopf,
        -1,  # dummy index
        vcat(x0, t0[1]),
        J_sc[1:N, 1:N+1],
        WeakRef(),  # dummy
    )
    param = set(solver.prob.param_axis2, de_prob.p, u_sc[end])

    return LimitCycleProblem(
        hopf_point,
        solver.opts;
        de_prob = remake(de_prob; p = param),
        param_axis = solver.prob.param_axis1,
        t_domain = (solver.prob.t_domain[1][1],
                    solver.prob.t_domain[2][1]),
        # phase_space = solver.prob.phase_space,
        kwargs...)
end
