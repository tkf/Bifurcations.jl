using ..ArrayUtils: _eigvals, _eig, canonicalize
using ..FDiffUtils: vderiv2, deriv2, deriv3
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

    if cache.prev_lyapunov_coefficient * cache.lyapunov_coefficient < 0
        return PointTypes.bautin
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

sn_quadratic_coefficient(cache) =
    sn_quadratic_coefficient(timekind(cache), cache)

function sn_quadratic_coefficient(tkind, wrapper)
    cache = as(wrapper, Cachish)
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

function first_lyapunov_coefficient(wrapper)
    cache = as(wrapper, Cachish)
    x0 = ds_state(cache)
    A = ds_jacobian(cache)

    # TODO: find a nicer way to compute left/right eigenvectors
    rvals, rvecs = _eig(A)
    _, ri = findmin(abs.(real.(rvals)))
    ω₀ = imag(rvals[ri])
    q = rvecs[:, ri]
    if ω₀ < 0
        ω₀ *= -1
        q = conj(q)
    end
    lvals, lvecs = _eig(A')
    _, li = findmin(abs.(lvals .+ im * ω₀))
    p = lvecs[:, li]

    q = cast_container(typeof(x0), q)
    p = cast_container(typeof(x0), p)
    q = canonicalize(q)     # s.t. q ⋅ q ≈ 1 and real(q) ⋅ imag(q) ≈ 0
    p /= p ⋅ q              # s.t. p ⋅ q ≈ 1
    qR = real(q)
    qI = imag(q)
    pR = real(p)
    pI = imag(q)

    a = vderiv2(t -> ds_f((@. t * qR + x0), cache))
    b = vderiv2(t -> ds_f((@. t * qI + x0), cache))
    c = 1/4 * vderiv2(t -> (ds_f((@. t * (qR + qI) + x0), cache) -
                            ds_f((@. t * (qR - qI) + x0), cache)))
    r = A \ (@. a + b)  # = A⁻¹ B(q,q̅)
    s = (@. 2im * ω₀ - A) \ (@. a - b + 2im * c)  # = (2iω - A)⁻¹ B(q,q)
    sR = real(s)
    sI = imag(s)

    σ₁ = 1/4 * deriv2(t -> (pR ⋅ ds_f((@. t * (qR + r) + x0), cache) -
                            pR ⋅ ds_f((@. t * (qR - r) + x0), cache)))
    σ₂ = 1/4 * deriv2(t -> (pI ⋅ ds_f((@. t * (qR + r) + x0), cache) -
                            pI ⋅ ds_f((@. t * (qR - r) + x0), cache)))
    Σ₀ = σ₁ + σ₂

    δ₁ = 1/4 * deriv2(t -> (pR ⋅ ds_f((@. t * (qR + sR) + x0), cache) -
                            pR ⋅ ds_f((@. t * (qR - sR) + x0), cache)))
    δ₂ = 1/4 * deriv2(t -> (pR ⋅ ds_f((@. t * (qR + sI) + x0), cache) -
                            pR ⋅ ds_f((@. t * (qR - sI) + x0), cache)))
    δ₃ = 1/4 * deriv2(t -> (pI ⋅ ds_f((@. t * (qR + sI) + x0), cache) -
                            pI ⋅ ds_f((@. t * (qR - sI) + x0), cache)))
    δ₄ = 1/4 * deriv2(t -> (pI ⋅ ds_f((@. t * (qR + sR) + x0), cache) -
                            pI ⋅ ds_f((@. t * (qR - sR) + x0), cache)))
    Δ₀ = δ₁ + δ₂ + δ₃ - δ₄

    γ₁ = deriv3(t -> pR ⋅ ds_f((@. t * qR + x0), cache))
    γ₂ = deriv3(t -> pI ⋅ ds_f((@. t * qI + x0), cache))
    γ₃ = deriv3(t -> (@. pR + pI) ⋅ ds_f((@. t * (qR + qI) + x0), cache))
    γ₄ = deriv3(t -> (@. pR - pI) ⋅ ds_f((@. t * (qR - qI) + x0), cache))
    Γ₀ = (2/3) * (γ₁ + γ₂) + (1/6) * (γ₃ + γ₄)  # = Re[<p, C(q,q,q̅)>]

    return (Γ₀ - 2Σ₀ + Δ₀) / 2ω₀
end

function testfn(pvtype::Val{PointTypes.cusp},
                tkind,
                ::SaddleNodeCont,
                prob_cache, u, J, L, Q)
    return sn_quadratic_coefficient(tkind, FakeCache(prob_cache, u, J))
end

function testfn(pvtype::Val{PointTypes.bautin},
                ::Continuous,
                ::HopfCont,
                prob_cache, u, J, L, Q)
    return first_lyapunov_coefficient(FakeCache(prob_cache, u, J))
end

function testfn(pvtype::Val{PointTypes.bogdanov_takens},
                ::Continuous,
                ::HopfCont,
                prob_cache, u, J, L, Q)
    return u[end - 2]  # = w = angular velocity
end

function testfn(::Union{Val{PointTypes.bogdanov_takens},
                        Val{PointTypes.fold_hopf}},
                ::Continuous, ckind::SaddleNodeCont,
                prob_cache, u, J, L, Q)
    vals = sort_by_abs_real!(_eigvals(ds_jacobian(ckind, J)))
    return real(vals[2])
end

function testfn(::Union{Val{PointTypes.hopf_hopf},
                        Val{PointTypes.fold_hopf}},
                ::Continuous, ckind::HopfCont,
                prob_cache, u, J, L, Q)
    return det(ds_jacobian(ckind, J))
end
