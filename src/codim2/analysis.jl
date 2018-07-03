using ..ArrayUtils: _eigvals, _eig
using ..Codim1: ds_eigvals
import ..Codim1: testfn

function sort_by_abs_real!(vals)
    sort!(vals, by=abs ∘ real)
    return vals
end

function guess_point_type(::Any, ::Discrete, cache, opts)
    return PointTypes.none
end

function guess_point_type(::SaddleNodeCont, ::Continuous, cache, opts)
    if cache.prev_quadratic_coefficient * cache.quadratic_coefficient < 0
        return PointTypes.cusp
    end

    if length(cache.eigvals) >= 2
        ev0 = sort_by_abs_real!(copy(cache.eigvals))[2]
        ev1 = sort_by_abs_real!(copy(cache.prev_eigvals))[2]
        if real(ev0) * real(ev1) <= 0
            if mean(imag.((ev0, ev1))) <= opts.atol
                return PointTypes.bogdanov_takens
            else
                return PointTypes.fold_hopf
            end
        end
    end

    return PointTypes.none
end

function guess_point_type(::HopfCont, ::Continuous, cache, opts)
    w = as(cache, ContinuationCache).u[end - 2]
    if w <= 0
        # Relying on that it was started positive.
        # See: [[./problem.jl::w0 < 0]]
        return PointTypes.bogdanov_takens
    end

    if length(cache.eigvals) >= 3
        # TODO: directly compute determinant and avoid computing
        # eigenvalues all the time
        d0 = prod(cache.prev_eigvals)
        d1 = prod(cache.eigvals)
        if d0 * d1 <= 0
            ev0 = sort_by_abs_real!(copy(cache.eigvals))[3]
            ev1 = sort_by_abs_real!(copy(cache.prev_eigvals))[3]
            if mean(imag.((ev0, ev1))) <= opts.atol
                return PointTypes.fold_hopf
            else
                return PointTypes.hopf_hopf
            end
        end
    end

    return PointTypes.none
end

function right_eigvec(::Discrete, J)
    vals, vecs = _eig(J)
    val, idx = findmax(abs.(vals))
    # @assert val ≈ 1
    v0 = vecs[:, idx]
    return v0 ./ norm(v0)
end

function right_eigvec(::Continuous, J)
    vals, vecs = _eig(J)
    val, idx = findmax(real.(vals))
    # @assert val ≈ 0
    v0 = vecs[:, idx]
    return v0 ./ norm(v0)
end
# TOOD: use eigs (depending on size(J)?)

left_eigvec(tkind, J) = right_eigvec(tkind, J')  # TODO: optimize

function sn_quadratic_coefficient(tkind, ::ImmutableState, wrapper)
    cache = as(wrapper, ContinuationCache)
    x0 = ds_state(cache)
    J = ds_jacobian(cache)
    q = right_eigvec(tkind, J)
    p = left_eigvec(tkind, J)
    q = cast_container(typeof(x0), q)
    p = cast_container(typeof(x0), p)
    second_derivative = ForwardDiff.hessian(
        t -> p ⋅ ds_f((@. t[1] * q + x0), cache),
        SVector(zero(eltype(x0))),
        # TODO: store config
    )[1, 1]
    return second_derivative / 2
end

# I need to access the cache for this:
#=
function testfn(pvtype::Val{PointTypes.cusp},
                tkind,
                ckind::SaddleNodeCont,
                u, J, L, Q)
end
=#

function testfn(pvtype::Val{PointTypes.bogdanov_takens},
                ::Continuous,
                ::HopfCont,
                u, J, L, Q)
    return u[end - 2]  # = w = angular velocity
end

function testfn(::Union{Val{PointTypes.bogdanov_takens},
                        Val{PointTypes.fold_hopf}},
                ::Continuous, ckind::SaddleNodeCont,
                u, J, L, Q)
    vals = sort_by_abs_real!(_eigvals(ds_jacobian(ckind, J)))
    return real(vals[2])
end

function testfn(::Union{Val{PointTypes.hopf_hopf},
                        Val{PointTypes.fold_hopf}},
                ::Continuous, ckind::HopfCont,
                u, J, L, Q)
    return det(ds_jacobian(ckind, J))
end
