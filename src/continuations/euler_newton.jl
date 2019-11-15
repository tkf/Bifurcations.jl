using Parameters: @with_kw

# TODO: Cleanup!  The code is very ugly at the moment because I
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
                                 uType, HType, JType, QType, hType,
                                 } <: AbstractContinuationCache{PC}
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
        _similar(u0, N, N),      # Q
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
    h0::Float64 = 0.01
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
    nominal_angle_rad::Float64 = 2π * (10 / 360)
    start_from_nearest_root::Bool = false
    bidirectional_first_sweep::Bool = true
    verbose::Bool = false
end

rawtangent(Q) = vec(rawtangentmat(Q))
function rawtangentmat(Q)
    if Q isa StaticArray
        return bottomrow(Q)
    elseif Q isa Adjoint  # CuArray takes this path
        x′ = _similar(Q, size(Q, 1), 1)
        fill!(x′, false)
        x′[end, 1] = 1
        # return (Q' * x′)'  # `vec(x::Adjoint{_, <:CuArray})` does not work
        return conj(reshape(Q' * x′, 1, :))
    else
        x = _similar(Q, 1, size(Q, 1))
        fill!(x, false)
        # x .= (x .* false .+ CartesianIndices(x)) .== Ref(size(x))  # InvalidIRError
        x[1, end] = 1
        rmul!(x, Q)
        return x
    end
end

function tangent(L, Q)
    tJ = rawtangent(Q)
    if _det(Q) * _det(popbottomright(L)) < 0
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
    L, Q = _lq!(A)
    return tangent(L, Q)
end

function corrector_step!(H::HType,
                         J::JType,
                         Q::QType,
                         v::vType,
                         prob_cache) where {
                             HType <: AbstractVector,
                             JType <: AbstractMatrix,
                             QType <: AbstractMatrix,
                             vType <: AbstractVector,
                         }
    H, J = residual_jacobian!(H, J, v, prob_cache)
    A = vcat(J, _zeros(J, 1, size(J, 2)))  # TODO: improve
    L, Q = _lq!(A)

    # Following block is a workaround to make test_predator_prey.jl
    # etc. work.  It was working accidentally since an equivalent of
    # `popbottomright(L) \ H` below was calling non-StaticArrays
    # method.  (Although this may not be specific to StaticArrays; it
    # seems test_vs_svector.jl need this, too.)
    if isalmostzero(H, eps(eltype(H))) && abs(det(L)) < eps(eltype(L))
        # Do nothing when there is no need for correction (and
        # trying to do so throws).
        w = v
        dv = zero(v)
        return (w :: vType,
                dv,
                H :: HType,
                L,
                Q,
                J :: JType)
    end

    y = vcat(popbottomright(L) \ H, _zeros(J, 1))
    dv = Q' * y
    w = v - dv
    return (w :: vType,
            dv,
            H :: HType,
            L,
            Q,
            J :: JType)
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

    local v, L
    for _ in 1:opts.max_adaptations
        # predictor
        v = u .+ direction * h .* tJ

        # corrector
        v, dv, H, L, Q, J = corrector_step!(H, J, Q, v, prob_cache)
        n1 = norm(dv)
        tJv = tangent(L, Q)
        angle = acos(min(abs(tJ ⋅ tJv), 1))  # TODO: should I use min?
        @debug "maximum(abs, H) = $(maximum(abs, H))"

        if all(abs.(H) .< min(2 * eps(eltype(H)), atol))
            @debug "corrector_skipped: The first correction is too close to the zero point."
            # The first correction is too close to the zero point.
            # Then the fractions for step adaptation would not be
            # possible to reliably calculated so let's skip them.
            # Also, to get out of this "tricky" region faster, let's
            # increase `h` slightly.
            h = h * 2  # TODO: don't hard code
            @goto corrector_skipped
        end

        v, dv, H, L, Q, J = corrector_step!(H, J, Q, v, prob_cache)
        n2 = norm(dv)
        @debug "maximum(abs, H) = $(maximum(abs, H))"

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

        @debug(
            "Step adaptation",
            # Print some useful statistics:
            n2 / n1,
            n1,
            angle,
            f_contraction,
            f_distance,
            f_angle,
            f0
        )

        h = h / f
        if h < opts.h_min
            break
        elseif isalmostzero(H, atol)
            @debug "corrector_skipped: isalmostzero(H, atol)"
            @goto corrector_skipped
        elseif f0 <= 2
            @debug "adaptation_success"
            cache.adaptation_success = true
            @goto adaptation_success
        end
    end
    # TODO: redesign
    error("""Failed to adapt steplength h.
        h = $h
        h < h_min = $(h < opts.h_min)
        max_adaptations = $(opts.max_adaptations)
        """)
    # step adaptation failed
    return

    @label adaptation_success
    # corrector (again)
    for _ in 3:opts.max_corrector_steps
        if isalmostzero(H, atol)
            cache.corrector_success = true
            @goto corrector_success
        end
        v, _, H, L, Q, J = corrector_step!(H, J, Q, v, prob_cache)
    end
    error("Failed in corrector loop.")  # TODO: redesign
    # corrector failed
    return

    # TODO: consider if I can avoid this special case
    @label corrector_skipped
    cache.adaptation_success = true
    cache.corrector_success = true

    @label corrector_success
    cache.u = v
    cache.J = J
    cache.h = h

    tJv = tangent(L, Q)
    if tJ ⋅ tJv < 0
        @debug "switching direction: $(cache.direction) -> $(-cache.direction)"
        cache.direction *= -1
        cache.simple_bifurcation = true
    end
    return v, h
end
# TODO: make it possible to compile away `@debug`s?


function nearest_root!(cache::ContinuationCache,
                       opts::ContinuationOptions,
                       v = cache.u)
    prob_cache = cache.prob_cache
    H = cache.H
    J = cache.J
    Q = cache.Q
    rtol = opts.rtol
    atol = opts.atol

    cache.corrector_success = false
    for _ in 1:opts.max_corrector_steps
        v, _, H = corrector_step!(H, J, Q, v, prob_cache)
        if isalmostzero(H, atol)
            @goto corrector_success
        end
    end
    error("Failed in corrector loop.")  # TODO: redesign

    @label corrector_success
    cache.corrector_success = true
    return v
end
