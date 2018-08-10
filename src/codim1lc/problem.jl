using ForwardDiff
using Parameters: @with_kw, @unpack
using Setfield: set
using StaticArrays: SVector, SMatrix, Size

using Jacobi: lagrange, zgj, wgj
using ..PolyUtils: dlagrange

import ..Continuations: get_prob_cache, get_u0, residual!, residual_jacobian!,
    isindomain
using ..BifurcationsBase: MutableState, TimeKind

abstract type
    Codim1LCProblem{skind, tkind} <: BifurcationProblem{skind, tkind}
end

struct LimitCycleProblem{
        skind <: MutableState,  # <: StateKind,
        tkind <: TimeKind,
        XS <: AbstractMatrix,
        X,
        R <: Real,
        P,
        L,
        } <: Codim1LCProblem{skind, tkind}
    statekind::skind
    timekind::tkind
    xs0::XS
    l0::R  # period
    t0::R  # parameter  (TODO: rename parameter to something else)
    period_bound::Tuple{R, R}
    t_domain::Tuple{R, R}
    phase_space::Tuple{X, X}
    diameter_bound::Tuple{R, R}
    time_offset::R
    de_prob::P
    param_axis::L
    # TODO: Maybe move `num_mesh` and `degree` to solver options:
    num_mesh::Int
    degree::Int
end

function LimitCycleProblem(;
        xs0,
        l0,
        t0,
        period_bound = (0.0, Inf),
        t_domain = (typemin(typeof(t0)), typemax(typeof(t0))),
        phase_space = (typemin(eltype(xs0)),
                       typemax(eltype(xs0))),
        diameter_bound = (1e-5, Inf),
        time_offset = 0.0,
        de_prob,
        param_axis,
        num_mesh,
        degree,
        )
    return LimitCycleProblem(
        MutableState(),
        timekind(de_prob),
        xs0,
        l0,
        t0,
        period_bound,
        t_domain,
        phase_space,
        diameter_bound,
        time_offset,
        de_prob,
        param_axis,
        num_mesh,
        degree,
    )
end

dim_state(prob::LimitCycleProblem) = size(prob.xs0, 1)
degree(prob::LimitCycleProblem) = prob.degree
num_mesh(prob::LimitCycleProblem) = prob.num_mesh

# ----------------------------------------------------------------------- cache

struct LimitCycleCache{P, VRef, MV, ML, VV} <: AbstractProblemCache{P}
    prob::P
    reference::VRef
    dv::MV
    lagrange_polynomial_vals::ML
    lagrange_polynomial_driv::ML
    gauss_quadrature_weight::VV
end

# Indices:
#   i  ∈  1:m    :  collocation/evaluation/Gauss point/node index; ζᵢ
#   k  ∈  1:m+1  :  "sample" node index; x(τₖ)
# Arrays:       Index:    Size:
#   lpv, lpd    [i, k]    (m, m + 1)
#   weight      [i]       (m,)
# Abbreviations:
#   lpv     :  lagrange_polynomial_vals
#   lpd     :  lagrange_polynomial_driv
#   weight  :  gauss_quadrature_weight

function LimitCycleCache(prob)
    n = dim_state(prob)
    m = degree(prob)
    reference = reshape(copy(prob.xs0), length(prob.xs0))
    dv = similar(prob.xs0, (n, m))

    # Compute τ: equi-distant points (within a mesh)
    interval = 1 / num_mesh(prob)
    dτ = interval / m
    τ = range(0, dτ, m + 1)
    @assert τ[end] ≈ interval

    ζ = zgj(m, 0, 0) .* interval  # Gauss nodes
    gauss_quadrature_weight = wgj(ζ, 0, 0)
    lagrange_polynomial_vals = [lagrange(k, ζᵢ, τ) for ζᵢ in ζ, k in 1:m+1]
    lagrange_polynomial_driv = [dlagrange(k, ζᵢ, τ) for ζᵢ in ζ, k in 1:m+1]

    return LimitCycleCache(
        prob,
        reference,
        SMatrix{n, m}(dv),
        SMatrix{m, m + 1}(lagrange_polynomial_vals),
        SMatrix{m, m + 1}(lagrange_polynomial_driv),
        SVector{m}(gauss_quadrature_weight),
    )
end
# TODO: dv is only used to hold the dimension of the dynamical system.
# Remove it.

Base.@pure dim_state(cache::LimitCycleCache) = Size(cache.dv)[1]
Base.@pure degree(cache::LimitCycleCache) = Size(cache.dv)[2]
num_mesh(cache::LimitCycleCache) = cache.prob.num_mesh

function set_reference!(cache::LimitCycleCache, u)
    copy!(cache.reference, @view u[1:length(cache.reference)])
end

# ------------------------------------------------------ continuation interface

get_prob_cache(prob::LimitCycleProblem) = LimitCycleCache(prob)

get_u0(prob::LimitCycleProblem) = vcat(
    view(prob.xs0, :),
    prob.l0,
    prob.t0,
)

u_idx_period(cache) = length(cache.prob.xs0) + 1
u_idx_param(cache::LimitCycleCache) = u_idx_period(cache) + 1
H_idx_phase_condition(cache) = length(cache.prob.xs0) + 1

isindomain(u, cache::LimitCycleCache) = isindomain_lc(u, cache)

function isindomain_lc(u, wrapper)
    cache = as(wrapper, LimitCycleCache)
    tmin, tmax = wrapper.prob.t_domain
    all(tmin .<= u[u_idx_param(wrapper)] .<= tmax) || return false

    period_min, period_max = cache.prob.period_bound
    (period_min <= u[u_idx_period(wrapper)] <= period_max) || return false

    xmin, xmax = cache.prob.phase_space
    n = dim_state(cache)
    for i in 1:degree(cache) * num_mesh(cache)
        x = @view u[n * (i - 1) + 1:n * i]
        all(xmin .<= x .<= xmax) || return false
    end

    xs = reshape((@view u[1:length(cache.prob.xs0)]), size(cache.prob.xs0))
    in_diameter_bound(xs, cache.prob.diameter_bound) || return false

    return true
end

function in_diameter_bound(xs::AbstractMatrix, diameter_bound)
    diameter_min, diameter_max = diameter_bound
    if diameter_max == typemax(diameter_max) && diameter_min == 0.0
        return true
    elseif diameter_max == typemax(diameter_max)
        return any(d -> d >= diameter_min, iter_diameter_candidates(xs))
    elseif diameter_min == 0.0
        return all(d -> d <= diameter_max, iter_diameter_candidates(xs))
    else
        d = maximum(iter_diameter_candidates(xs))
        return (diameter_min <= d <= diameter_max)
    end
end

iter_diameter_candidates(xs::AbstractMatrix) =
    (norm(@views xs[:, i] .- xs[:, j])
     for i in 1:size(xs, 2)
     for j in i + 1:size(xs, 2))

function diameter(u::AbstractVector, prob::LimitCycleProblem)
    xs = reshape((@view u[1:length(prob.xs0)]), size(prob.xs0))
    return maximum(iter_diameter_candidates(xs))
end

# ------------------------------------------------------------------- workspace

@inline function diff_range(cache, j, i)
    n = dim_state(cache)
    m = degree(cache)
    offset = n * m * (j - 1) + n * (i - 1)
    return offset + 1:offset + n
end

# TODO: get_samples --- better name than "samples"?
@generated function get_samples(
        cache::LimitCycleCache{P, VRef, MV, ML, VV},
        u, j,
        ) where {P, VRef, MV, ML, VV}
    n = Size(MV)[1]
    m = Size(MV)[2]

    idxs_rest = [:($(n * m) * (j - 1) + $i) for i in  1:(n * (m + 1))]
    idxs_last = [
        idxs_rest[1:end - n]
        1:n
    ]
    Mat = SMatrix{n, m + 1, eltype(u)}

    quote
        if j == num_mesh(cache)
            @inbounds return $Mat($([:(u[$i]) for i in idxs_last]...))
        else
            @inbounds return $Mat($([:(u[$i]) for i in idxs_rest]...))
        end
    end
end
# TODO: Stop hard-coding output type here. Assuming that m is always
# small (for StaticArrays to be effective) is probably fine but
# assuming n to be small is not OK.  But I need to optimize a lot for
# Bifurcations.jl to be useful for non-small n anyway.

# Indices:
#   i  ∈  1:m    :  collocation/evaluation/Gauss point/node index; ζᵢ
#   k  ∈  1:m+1  :  "sample" node index; x(τₖ)
#   p  ∈  1:n    :  phase space dimension
# Arrays:       Index:    Size:
#   x, dx       [p, i]    (n, m)
#   lpv, lpd    [i, k]    (m, m + 1)
#   xτ          [p, k]    (n, m + 1)

@inline function collocation(cache, u, j)
    lpv = cache.lagrange_polynomial_vals  # ℓᵀ; ℓᵢₖ = [ℓᵀ]ₖᵢ = ℓₖ(ζᵢ)
    lpd = cache.lagrange_polynomial_driv
    xτ = get_samples(cache, u, j) # xₚ(τₖ)
    x = xτ * lpv'                    # x(ζ) = x(τ) * ℓ
    dx = xτ * lpd'                   # x'(ζ) = x'(τ) * ℓ
    return x, dx
end
# TODO: store ℓ (not ℓᵀ) in lpv? (Why not?)

@inline function reference_diff(cache, j)
    dv = cache.dv
    lpd = cache.lagrange_polynomial_driv
    vτ = get_samples(cache, cache.reference, j)
    dv = vτ * lpd'  # x'(ζ) = x'(τ) * ℓ
    return dv
end

"""
    turns(num_mesh::Integer, degree::Integer, j::Integer)
    turns(cache, j::Integer)

Return an iterable of scaled times ([turn]) for `j`-th mesh.

[turn]: https://en.wikipedia.org/wiki/Turn_(geometry)
"""
turns(cache, j) = turns(num_mesh(cache), degree(cache), j)

function turns(num_mesh::Integer, degree::Integer, j::Integer)
    interval = 1 / num_mesh
    dt = interval / degree
    return range(interval * (j - 1), dt, degree)
end

# ------------------------------------------------------------------- residual!

function residual!(H, u, cache::LimitCycleCache)
    prob = cache.prob
    q = set(prob.param_axis, prob.de_prob.p, u[u_idx_param(cache)])
    return residual_lc!(H, u, q, cache)
end

function residual_lc!(H, u, q, cache::LimitCycleCache)
    prob = cache.prob

    l = u[u_idx_period(cache)]
    weight = cache.gauss_quadrature_weight

    phase_condition = zero(eltype(H))
    for j in 1:num_mesh(cache)
        x, dx = collocation(cache, u, j)
        dv = reference_diff(cache, j)
        @inbounds for (i, t) in enumerate(turns(cache, j))
            r = diff_range(cache, j, i)
            f = @view H[r]
            prob_time = t * l + prob.time_offset
            prob.de_prob.f(f, x[:, i], q, prob_time)

            @. H[r] = dx[:, i] - l * f
            phase_condition += weight[i] * (x[:, i] ⋅ dv[:, i])
        end
    end
    H[H_idx_phase_condition(cache)] = phase_condition

    return H
end
# Note: phase_condition above is integrated at collocation points
# (Gauss nodes) ζᵢ, not the Lagrange nodes τₖ, as mentioned in
# Kuznetsov (1998).

# ---------------------------------------------------------- residual_jacobian!

function residual_jacobian!(H, J, u, cache::LimitCycleCache)
    # TODO: implement residual_jacobian! manually
    ForwardDiff.jacobian!(
        J,
        (y, x) -> residual!(y, x, cache),
        H,  # y
        u,  # x
    )
    return (H, J)
end
