using Parameters: @with_kw
using StaticArrays: SMatrix, SVector

# TODO: Cleanup!  The code is the very ugly at the moment because I
# tried to support out-of-place and in-place algorithm in the same
# code (and failing to do so).  I should:
#
# * Add some tests for out-of-place code.
# * Add some benchmarks, especially for in-place code.
# * Reduce type-instability, especially for in-place code.
# * Find the right abstraction and then finally clean it up!

"""
Cache for Euler-Newton continuation method.

See [`AbstractContinuationProblem`](@ref) for the mathematical setup.

# Fields
- `prob_cache`
- `u` (size: `(N,)`)
- `H` (size: `(N - 1,)`) ``= H(u)``
- `J` (size: `(N - 1, N)`) ``= \\partial H / \\partial u``
- `Q` (size: `(N - 1, N)`): temporary array for the QR decomposition
- `h::Real`: step size
- `direction::Int`: +1 or -1
- `corrector_success::Bool`
- `adaptation_success::Bool`
- `simple_bifurcation::Bool`
"""
mutable struct ContinuationCache{PC <: AbstractProblemCache,
                                 uType, HType, JType, QType, hType}
    prob_cache::PC
    u::uType
    H::HType
    J::JType
    Q::QType
    h::hType
    direction::Int
    corrector_success::Bool
    adaptation_success::Bool
    simple_bifurcation::Bool
end

function ContinuationCache(prob_cache::AbstractProblemCache,
                           h::Real, direction::Int = 1)
    u0 = get_u0(prob_cache.prob)
    N = length(u0)
    return ContinuationCache(
        prob_cache,
        u0,
        _similar(u0, N - 1),     # H
        _similar(u0, N - 1, N),  # J
        _similar(u0, N - 1, N),  # Q
        h,
        direction,
        false,
        false,
        false,
    )
end

ContinuationCache(prob::AbstractContinuationProblem, args...) =
    ContinuationCache(get_prob_cache(prob), args...)


@with_kw struct ContinuationOptions
    direction::Int = 1
    h0::Float64 = 1.0
    h_min::Float64 = 1e-6
    h_zero::Float64 = 1e-6
    rtol::Float64 = 0.01
    atol::Float64 = 1e-6
    max_samples::Int = 100
    max_adaptations::Int = 100
    max_corrector_steps::Int = 100
    max_branches::Int = 10
    max_misc_steps::Int = 100   # TODO: remove
    nominal_contraction::Float64 = 0.8
    nominal_distance::Float64 = 0.1
    nominal_angle_rad::Float64 = 2π * (30 / 360)
end


function tangent(L, Q)
    tJ = Q[:, end]
    if det(Q) * det(@view L[1:end-1, 1:end-1]) < 0
        tJ *= -1
    end
    return tJ
end

function current_tangent(cache::ContinuationCache,
                         opts::ContinuationOptions)
    prob_cache = cache.prob_cache
    u = cache.u
    H = cache.H
    J = cache.J
    Q = cache.Q

    H, J = residual_jacobian!(H, J, u, prob_cache)
    A = vcat(J, _zeros(J, 1, size(J, 2)))  # TODO: improve
    L, Q = lq!(Q, A)
    return tangent(L, Q)
end

function corrector_step!(H, J, Q, v, prob_cache)
    H, J = residual_jacobian!(H, J, v, prob_cache)
    A = vcat(J, _zeros(J, 1, size(J, 2)))  # TODO: improve
    L, Q = lq!(Q, A)
    y = _A_ldiv_B!((@view L[1:end-1, 1:end-1]), H)
    dv = (@view Q[:, 1:end-1]) * y
    w = v - dv
    return (w, dv, H, L, Q, J)
end

function predictor_corrector_step!(cache::ContinuationCache,
                                   opts::ContinuationOptions)
    predictor_corrector_step!(cache, opts,
                              cache.u,
                              current_tangent(cache, opts),
                              cache.h,
                              cache.direction)
end

function predictor_corrector_step!(cache::ContinuationCache,
                                   opts::ContinuationOptions,
                                   u, tJ, h, direction)
    prob_cache = cache.prob_cache
    H = cache.H
    J = cache.J
    Q = cache.Q
    rtol = opts.rtol
    atol = opts.atol

    cache.corrector_success = false
    cache.adaptation_success = false
    cache.simple_bifurcation = false

    for _ in 1:opts.max_adaptations
        # predictor
        v = u .+ direction * h .* tJ

        # corrector
        v, dv, H, L, Q, _ = corrector_step!(H, J, Q, v, prob_cache)
        n1 = norm(dv)
        tJv = tangent(L, Q)
        angle = acos(min(abs(tJ ⋅ tJv), 1))  # TODO: should I use min?

        v, dv, H, L, Q, _ = corrector_step!(H, J, Q, v, prob_cache)
        n2 = norm(dv)

        # step adaptation
        f_contraction =
            sqrt(n2 / n1 / opts.nominal_contraction)   # √ κ(u,h) / κ̃
        f_distance = sqrt(n1 / opts.nominal_distance)  # √ δ(u,h) / δ̃
        f_angle = angle / opts.nominal_angle_rad       #   α(u,h) / α̃
        f0 = max(
            zero_if_nan(f_contraction),
            zero_if_nan(f_distance),
            zero_if_nan(f_angle),
        )
        f = max(min(f0, 2), 1/2)
        h = h / f
        if h < opts.h_min
            return
        # elseif isalmostzero(H, rtol, atol)
        #     # If close enough to the solution, let it pass?  Should I?
        elseif f0 > 2
            continue
        end

        # corrector (again)
        for _ in 3:opts.max_corrector_steps
            if isalmostzero(H, rtol, atol)
                cache.corrector_success = true
                break
            end
            v, _, H, L, Q, _ = corrector_step!(H, J, Q, v, prob_cache)
        end
        if ! cache.corrector_success
            return
        end

        cache.u = v
        cache.h = h
        cache.adaptation_success = true

        tJv = tangent(L, Q)
        if tJ ⋅ tJv < 0
            cache.direction *= -1
            cache.simple_bifurcation = true
        end
        return v, h
    end
end
