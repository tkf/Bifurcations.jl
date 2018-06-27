using ..Codim1: ds_eigvals

function guess_point_type(::Any, ::Discrete, cache, opts)
    return PointTypes.none
end

function guess_point_type(::SaddleNodeCont, ::Continuous, cache, opts)
    if cache.prev_quadratic_coefficient * cache.quadratic_coefficient < 0
        return PointTypes.cusp
    end

    sign_changed, is_real = eigval_changes(cache, opts, 2)
    if sign_changed
        if is_real
            return PointTypes.bogdanov_takens
        else
            return PointTypes.fold_hopf
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

    # i=3 to skip the conjugat; probably not robust...  # TODO: fix
    sign_changed, is_real = eigval_changes(cache, opts, 3)
    if sign_changed
        if is_real
            return PointTypes.fold_hopf
        else
            return PointTypes.hopf_hopf
        end
    end

    return PointTypes.none
end

function eigval_changes(cache, opts, i)
    if length(cache.eigvals) >= i
        ev0 = cache.prev_eigvals[i]
        ev1 = cache.eigvals[i]
        if real(ev0) * real(ev1) < 0
            img = mean(abs.(imag.((ev0, ev1))))
            return true, img < opts.atol
        end
    end
    return false, false
end

function right_eigvec(::Discrete, J)
    vals, vecs = eig(J)
    val, idx = findmax(abs.(vals))
    # @assert val ≈ 1
    v0 = vecs[:, idx]
    return v0 ./ norm(v0)
end

function right_eigvec(::Continuous, J)
    vals, vecs = eig(J)
    val, idx = findmax(real.(vals))
    # @assert val ≈ 0
    v0 = vecs[:, idx]
    return v0 ./ norm(v0)
end
# TOOD: use eigs (depending on size(J)?)

left_eigvec(tkind, J) = right_eigvec(tkind, J')  # TODO: optimize

function sn_quadratic_coefficient(tkind, ::ImmutableState, wrapper)
    cache = as(wrapper, ContinuationCache)
    J = ds_jacobian(cache)
    q = right_eigvec(tkind, J)
    p = left_eigvec(tkind, J)
    x0 = ds_state(cache)
    second_derivative = ForwardDiff.hessian(
        t -> p ⋅ ds_f((@. t[1] * q + x0), cache),
        SVector(zero(x0)),
        # TODO: store config
    )[1, 1]
    return second_derivative / 2
end
