using ..Continuations: _normalize!
using ..ArrayUtils: _eig

other(::HopfCont) = SaddleNodeCont()
other(::SaddleNodeCont) = HopfCont()

orthogonalize(v, u) = v .- u .* ((u ⋅ v) / (u ⋅ u))


function BifurcationProblem(point::AbstractSpecialPoint,
                            solver::Codim2Solver;
                            t_domain = solver.prob.t_domain,
                            param_axis1 = solver.prob.param_axis1,
                            param_axis2 = solver.prob.param_axis2,
                            kwargs...)

    cd2_prob = solver.prob :: DiffEqCodim2Problem
    de_prob = cd2_prob.de_prob :: AbstractODEProblem

    # Manually cast.  The fact that `resolved.u` is not of type xtype
    # indicates type instability somewhere.  Until it is fixed, let's
    # silently cast always.  After type stability is established, it
    # has to be changed to emit warning before cast.
    xtype = typeof(de_prob.u0)

    resolved = resolve_point(point, solver)
    N = length(de_prob.u0)
    x0 = xtype(resolved.u[1:N])
    t0 = SVector{2}(resolved.u[end-1:end])
    @assert param_axis1 == solver.prob.param_axis1
    @assert param_axis2 == solver.prob.param_axis2

    @assert point.point_type in (PointTypes.bogdanov_takens,
                                 PointTypes.fold_hopf)
    if (contkind(solver) :: ContinuationKind) isa SaddleNodeCont
        next_ckind = HopfCont()
    else
        next_ckind = SaddleNodeCont()
    end

    @assert timekind(point) isa Continuous

    v0, w0 = next_eig(next_ckind, resolved)
    @assert length(v0) == length(x0)
    if next_ckind isa HopfCont
        v0 = v0 .+ 0im
    else
        v0 = real.(v0)
    end
    v0 = cast_container(xtype, v0)

    # Sanity checks:
    if next_ckind isa SaddleNodeCont
        @assert eltype(v0) <: Real
    else
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
        augmented_system = preferred_augsys(next_ckind),
        kwargs...)
end

next_eig(next_ckind::ContinuationKind, point::SpecialPoint) =
    next_eig(next_ckind, Val(point.point_type), point)

# Switching to Hopf continuation via Bogdanov-Takens from saddle-node
function next_eig(::HopfCont, ::Val{PointTypes.bogdanov_takens},
                  point)
    prev_ckind = SaddleNodeCont()
    v_prev = ds_eigvec(prev_ckind, point.u)
    vals, vecs = _eig(ds_jacobian(prev_ckind, point.J))

    if eltype(vals) <: Complex
        _, i0 = findmin(abs.(real.(vals)))
        v0 = vecs[:, i0]
        w0 = imag(vals[i0])
        if w0 < 0
            w0 = -w0
            v0 = conj(v0)
        end
        return canonicalize(v0), w0
    end

    i1, i2 = sortperm(vals, by=abs)
    v1 = vecs[:, i1]
    v2 = vecs[:, i2]
    if abs(v1 ⋅ v_prev) > abs(v2 ⋅ v_prev)
        v_img = v1
    else
        v_img = v2
    end
    v_next = canonicalize(@. v_prev + v_img * im)
    return v_next, zero(real(eltype(v_prev)))
end

# Switching to saddle-node continuation via Bogdanov-Takens from Hopf
function next_eig(::SaddleNodeCont, ::Val{PointTypes.bogdanov_takens},
                  point)
    prev_ckind = HopfCont()
    v_prev = ds_eigvec(prev_ckind, point.u)
    v_next = _normalize!(real.(v_prev))
    return v_next, zero(eltype(v_next))
end

# Switching to Hopf (saddle-node) continuation via Fold-Hopf from
# saddle-node (Hopf)
function next_eig(next_ckind::ContinuationKind,
                  ::Val{PointTypes.fold_hopf},
                  point)
    prev_ckind = other(next_ckind)
    vals, vecs = _eig(ds_jacobian(prev_ckind, point.J))
    i1, i2, i3 = sortperm(vals, by=abs∘real)
    ir, _, ii = sort((i1, i2, i3), by=i -> abs(imag(vals[i])))
    if next_ckind isa HopfCont
        return _normalize!(vecs[:, ii]), imag(vals[ii])
    else
        return _normalize!(vecs[:, ir]), imag(vals[ir])
    end
end
