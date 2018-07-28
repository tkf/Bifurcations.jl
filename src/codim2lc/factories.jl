using DiffEqBase: AbstractODEProblem, remake
using StaticArrays: SVector
using Setfield: get, set

using ..Continuations: init, solve!, residual_jacobian, ContinuationSolver,
    ZeroNotFoundError, FindZeroInputError
using ..BifurcationsBase: AbstractSpecialPoint, SpecialPoint, special_points
using ..Codim1
using ..Codim1: resolving_points
using ..Codim2
using ..Codim2: Codim2Solver, DiffEqCodim2Problem, first_lyapunov_coefficient
using ..Codim1LimitCycle: Codim1LCSolver, diameter

# See also:
# * [[../codim2/problem.jl::function BifurcationProblem]]
# * [[../codim2/switch.jl::function BifurcationProblem]]
# * [[../codim1lc/factories.jl::function LimitCycleProblem]]
function FoldLimitCycleProblem(
        point::AbstractSpecialPoint,
        solver::Codim2Solver;
        lc_solver_opts = [],
        kwargs...)
    @assert point.point_type == Codim2.PointTypes.bautin
    @assert timekind(point) isa Continuous

    solver.prob :: DiffEqCodim2Problem
    de_prob = solver.prob.de_prob :: AbstractODEProblem

    lc_prob = bautin_to_subcritical_lc(point, solver; kwargs...)
    # TODO: Stop solving sub-problems in the constructor:
    u_init = nearest_fold_lc(lc_prob; lc_solver_opts...)

    prob_cache = get_prob_cache(lc_prob)
    _H, J_init = residual_jacobian(u_init, prob_cache)

    n = length(lc_prob.xs0) + 1  # +1 to include period
    vals, vecs = eig(@view J_init[1:n, 1:n])
    _, idx0 = findmin(abs.(vals))  # closest to zero
    opts = as(solver, ContinuationSolver).opts
    if any(x -> abs(imag(x)) > opts.atol, @view vecs[:, idx0])
        @warn "Nearest-to-null vector is complex (eigval: $(vals[idx0]))"
    end
    v0 = normalize(real(@view vecs[:, idx0]))

    xs0 = reshape(u_init[1:length(lc_prob.xs0)], size(lc_prob.xs0))
    l0 = u_init[length(lc_prob.xs0) + 1]
    vs0 = reshape(v0[1:end-1], size(lc_prob.xs0))
    dl0 = v0[end]
    t0 = SVector(u_init[end],
                 get(solver.prob.param_axis2, lc_prob.de_prob.p))

    return FoldLimitCycleProblem(
        lc_prob;
        xs0 = xs0,
        l0 = l0,
        vs0 = vs0,
        dl0 = dl0,
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


function nearest_fold_lc(lc_prob; lc_solver_opts...)
    lc_solver = init(
        lc_prob;
        start_from_nearest_root = true,
        max_branches = 0,
        nominal_angle_rad = 2π * (1 / 360),  # TODO: should it be the default?
        lc_solver_opts...
    )
    solve!(lc_solver)

    PT = Codim1LimitCycle.PointTypes.saddle_node
    point_candidates = special_points(lc_solver, PT)
    if isempty(point_candidates)
        error("Cannot find fold bifurcation of the limit cycle.")
    end
    sw1_points = special_points(lc_solver.sol.sweeps[1], PT)
    sw2_points = special_points(lc_solver.sol.sweeps[2], PT)
    sweep_index = -1
    if length(sw2_points) == 0
        @assert length(sw1_points) > 0
        flc_point, = sw1_points
        sweep_index = 1
    elseif length(sw1_points) == 0
        @assert length(sw2_points) > 0
        flc_point, = sw2_points
        sweep_index = 2
    else
        # Choose increasing sweep
        cont_sweeps = as(lc_solver, ContinuationSolver).sol.sweeps
        at_most = (n, u) -> u[1:max(length(u), n)]
        ds1 = [diameter(u, lc_prob) for u in at_most(10, cont_sweeps[1].u)]
        ds2 = [diameter(u, lc_prob) for u in at_most(10, cont_sweeps[2].u)]
        mean_dd1 = mean(diff(ds1))
        mean_dd2 = mean(diff(ds2))
        if mean_dd1 < 0 && mean_dd2 < 0
            @warn "Both sweeps have decreasing diameters."
        end
        if mean_dd1 > mean_dd2
            flc_point, = sw1_points
            sweep_index = 1
        else
            flc_point, = sw2_points
            sweep_index = 2
        end
    end

    # TODO?: resolve flc_point
    u_init = (flc_point.u0 .+ flc_point.u1) ./ 2

    if length(point_candidates) > 1
        @warn """
        $(length(point_candidates)) fold bifurcations are found.
        Using the first point found in sweep $sweep_index:
            $(lc_prob.param_axis) : $(u_init[end])
            diameter : $(diameter(u_init, lc_prob))
        """
    end

    return u_init
end


function bautin_to_subcritical_lc(
        point::AbstractSpecialPoint,
        solver::Codim2Solver;
        kwargs...)

    @assert point.point_type == Codim2.PointTypes.bautin
    @assert timekind(point) isa Continuous

    solver.prob :: DiffEqCodim2Problem
    de_prob = solver.prob.de_prob :: AbstractODEProblem
    prob_cache = as(solver.cache, ContinuationCache).prob_cache

    if point isa SpecialPoint
        point = let
            if point.sweep.value == nothing
                error("""
                `point` is not associated with a sweep.
                point = $(point)
                 """)
            end
            # See: [[../base/solver.jl::SpecialPointInterval]]
            sweep = as(point.sweep.value, ContinuationSweep)
            u0 = sweep.u[point.point_index - 1]
            u1 = sweep.u[point.point_index]
            _H, J1 = residual_jacobian(u1, prob_cache)
            SpecialPointInterval(
                timekind(point),
                point.point_type,
                point.point_index,
                u0,
                u1,
                J1,
                WeakRef(sweep),
            )
        end
    else
        point :: SpecialPointInterval
    end

    # TODO: Specify target l1 (> 0) value and step into it using zero
    # finder (so that hopefully it's not too close to/far away from
    # the Bautin point).
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
