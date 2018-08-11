using DiffEqBase: remake, AbstractODEProblem
using ...Bifurcations: FixedPointBifurcationProblem, DiffEqWrapper
using ..ArrayUtils: canonicalize
using ..Codim1
using ..Codim1: Codim1Solver, AbstractSpecialPoint, resolve_point

abstract type Codim2Problem{skind, tkind} <: BifurcationProblem{skind, tkind} end

"""
    BifurcationProblem(point::AbstractSpecialPoint,
                       solver::Codim1Solver,
                       param_axis2::Lens,
                       t2_domain::Tuple)

Construct codimension-2 bifurcation problem given a bifurcation
`point`.
"""
function BifurcationProblem(point::AbstractSpecialPoint,
                            solver::Codim1Solver,
                            param_axis2::Lens,
                            t2_domain::Tuple;
                            kwargs...)

    cd1_prob = solver.prob :: FixedPointBifurcationProblem
    de_wrapper = cd1_prob.p :: DiffEqWrapper
    de_prob = de_wrapper.de_prob :: AbstractODEProblem

    param_axis1 = de_wrapper.param_axis
    t_domain = (
        (cd1_prob.t_domain[1], t2_domain[1]),
        (cd1_prob.t_domain[2], t2_domain[2]),
    )

    # Manually cast.  The fact that `resolved.u` is not of type xtype
    # indicates type instability somewhere.  Until it is fixed, let's
    # silently cast always.  After type stability is established, it
    # has to be changed to emit warning before cast.
    xtype = typeof(de_prob.u0)

    resolved = resolve_point(point, solver)
    x0 = xtype(resolved.u[1:end - 1])
    t0 = SVector(resolved.u[end], get(param_axis2, de_prob.p))

    vals, vecs = eigen(Codim1.ds_jacobian(resolved.J))
    # TOOD: use eigs (depending on size(J)?)
    if timekind(point) isa Discrete
        val, idx = findmax(@. abs(abs(vals) - 1))
        # @assert val ≈ 1
    else
        val, idx = findmin(@. abs(real(vals)))
        # @assert val ≈ 0
    end
    v0 = cast_container(xtype, vecs[:, idx])
    v0 = canonicalize(v0)
    w0 = imag(vals[idx])
    if eltype(v0) <: Complex && w0 < 0
        w0 = -w0
        v0 = conj(v0)
    end

    # Sanity checks:
    if point.point_type == Codim1.PointTypes.saddle_node
        @assert eltype(v0) <: Real
    elseif point.point_type == Codim1.PointTypes.hopf
        @assert eltype(v0) <: Complex
    end
    if xtype <: SVector
        v0 :: SVector
    end

    return DiffEqCodim2Problem(
        de_prob,
        param_axis1,
        param_axis2,
        t_domain;
        x0 = x0,
        v0 = v0,
        w0 = w0,
        t0 = t0,
        augmented_system = preferred_augsys(point),
        kwargs...)
end
