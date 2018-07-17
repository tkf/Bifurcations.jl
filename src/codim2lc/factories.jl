using DiffEqBase: AbstractODEProblem
using StaticArrays: SVector

using ..BifurcationsBase: AbstractSpecialPoint, SpecialPoint
using ..Codim1
using ..Codim1: resolve_point
using ..Codim2
using ..Codim2: Codim2Solver, DiffEqCodim2Problem

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

    # [[../codim2/switch.jl::Manually cast]]
    xtype = typeof(de_prob.u0)

    resolved = resolve_point(point, solver)
    N = length(de_prob.u0)
    x0 = xtype(resolved.u[1:N])
    t0 = SVector{2}(resolved.u[end-1:end])
    param_axis1 = solver.prob.param_axis1
    param_axis2 = solver.prob.param_axis2

    hopf_point = SpecialPoint(
        timekind(point),
        Codim1.PointTypes.hopf,
        -1,  # dummy index
        vcat(x0, t0[1]),
        point.J[1:N, 1:N+1],
        WeakRef(),  # dummy
    )
    lc_prob = LimitCycleProblem(
        hopf_point,
        solver.opts;
        de_prob = de_prob,
        param_axis = param_axis1,
        kwargs...)
    # [[../codim1lc/factories.jl::function LimitCycleProblem]]

    center = hopf_point.u[1:end - 1]
    vs0 = lc_prob.xs0 .- center
    normalize!(view(vs0, :))

    return FoldLimitCycleProblem(
        lc_prob;
        vs0 = vs0,
        dl0 = 0.0,  # ???
        t0 = t0,
        param_axis1 = param_axis1,
        param_axis2 = param_axis2,
        t_domain = solver.prob.t_domain,
    )
end

# I need to implement fold bifurcation detection first.
#=
function BifurcationProblem(point::AbstractSpecialPoint,
                            solver::Codim1LCSolver,
                            param_axis2::Lens,
                            t2_domain::Tuple;
                            )

    cd1_prob = solver.prob :: LimitCycleProblem
    de_prob = cd1_prob.de_prob :: AbstractODEProblem

    xs0 = reshape(point.u[end:end-2], size(prob.xs0))

    # Initialize perturbation direction
    vs0 = similar(xs0)
    error("TODO: implement initial perturbation direction vs0")

    l0 = point.u[end - 1]
    t0 = SVector(point.u[end], get(param_axis2, de_prob.p))
    t_domain = (
        SVector(cd1_prob.t_domain[1], t2_domain[1]),
        SVector(cd1_prob.t_domain[2], t2_domain[2]),
    )
    param_axis1 = cd1_prob.param_axis

    return FoldLimitCycleProblem(
        cd1_prob.statekind,
        cd1_prob.timekind,
        cd1_prob,
        xs0,
        vs0,
        l0,
        t0,
        t_domain,
        param_axis1,
        param_axis2,
    )
end
=#
