using ..Codim1: ds_eigvals

function guess_point_type(::Discrete, cache, quadratic_coefficient, opts)
    return PointTypes.none
end

function guess_point_type(::Continuous, cache, quadratic_coefficient, opts)
    if quadratic_coefficient * cache.quadratic_coefficient < 0
        return PointTypes.cusp
    end
    return PointTypes.none
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
    x0 = ds_state(cache.u)
    second_derivative = ForwardDiff.hessian(
        t -> p ⋅ ds_f((@. t[1] * q + x0), cache),
        SVector(zero(x0)),
        # TODO: store config
    )[1, 1]
    return second_derivative / 2
end
