using ForwardDiff
using Parameters: @with_kw, @unpack
using Setfield: set

using Jacobi: lagrange, zgj, wgj
using ..PolyUtils: dlagrange

using ..CompatUtils: @required
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
    time_offset::R
    de_prob::P
    param_axis::L
    # TODO: Maybe move `num_mesh` and `degree` to solver options:
    num_mesh::Int
    degree::Int
end

function LimitCycleProblem(;
        (@required xs0),
        (@required l0),
        (@required t0),
        period_bound = (0.0, Inf),
        t_domain = (typemin(typeof(t0)), typemax(typeof(t0))),
        phase_space = (typemin(eltype(xs0)),
                       typemax(eltype(xs0))),
        time_offset = 0.0,
        (@required de_prob),
        (@required param_axis),
        (@required num_mesh),
        (@required degree),
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

struct LimitCycleCache{P, V, M} <: AbstractProblemCache{P}
    prob::P
    reference::V
    dv::M
    lagrange_polynomial_vals::M
    lagrange_polynomial_driv::M
    gauss_quadrature_weight::V
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
        dv,
        lagrange_polynomial_vals,
        lagrange_polynomial_driv,
        gauss_quadrature_weight,
    )
end

dim_state(cache::LimitCycleCache) = size(cache.dv, 1)
degree(cache::LimitCycleCache) = size(cache.dv, 2)
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
    return true
end

# ------------------------------------------------------------------- workspace

struct LimitCycleWorkspace{C, M}
    cache::C
    x::M
    dx::M
end

function make_workspace(cache, u)
    x = similar(u, size(cache.dv))
    dx = similar(x)
    return LimitCycleWorkspace(cache, x, dx)
end

function diff_range(cache, j, i)
    n = dim_state(cache)
    m = degree(cache)
    offset = n * m * (j - 1) + n * (i - 1)
    return offset + 1:offset + n
end

function get_samples(cache, u, j)   # TODO: better name than "samples"?
    n = dim_state(cache)
    m = degree(cache)
    is_last = j == num_mesh(cache)

    if is_last
        dims = (n, m)
    else
        dims = (n, m + 1)
    end
    s = n * m * (j - 1) + 1
    e = s - 1 + prod(dims)
    xτ = reshape(view(u, s:e), dims)
    xτ0 = view(u, 1:n)

    return xτ, xτ0, is_last
end

# Indices:
#   i  ∈  1:m    :  collocation/evaluation/Gauss point/node index; ζᵢ
#   k  ∈  1:m+1  :  "sample" node index; x(τₖ)
#   p  ∈  1:n    :  phase space dimension
# Arrays:       Index:    Size:
#   x, dx       [p, i]    (n, m)
#   lpv, lpd    [i, k]    (m, m + 1)
#   xτ          [p, k]    (n, m + 1) or (n, m)
#   xτ0         [p]       (n,)

function collocation!(ws, u, j)
    @unpack x, dx = ws
    lpv = ws.cache.lagrange_polynomial_vals  # ℓᵀ; ℓᵢₖ = [ℓᵀ]ₖᵢ = ℓₖ(ζᵢ)
    lpd = ws.cache.lagrange_polynomial_driv
    xτ, xτ0, is_last = get_samples(ws.cache, u, j) # xₚ(τₖ)
    if ! is_last
        A_mul_Bt!(x, xτ, lpv)  # x(ζ) = x(τ) * ℓ
        A_mul_Bt!(dx, xτ, lpd) # x'(ζ) = x'(τ) * ℓ
    else
        # It is the last mesh.  So, the last node has to come from the
        # first node of the first mesh.
        A_mul_Bt!(x, xτ, @view lpv[:, 1:end-1])
        A_mul_Bt!(dx, xτ, @view lpd[:, 1:end-1])
        x .+= xτ0 * (@view lpv[:, end])'
        dx .+= xτ0 * (@view lpd[:, end])'
        # TODO: optimize
    end
    return x, dx
end

function reference_diff!(cache, j)
    dv = cache.dv
    lpd = cache.lagrange_polynomial_driv
    vτ, vτ0, is_last = get_samples(cache, cache.reference, j)
    if ! is_last
        A_mul_Bt!(dv, vτ, lpd) # x'(ζ) = x'(τ) * ℓ
    else
        A_mul_Bt!(dv, vτ, @view lpd[:, 1:end-1])
        dv .+= vτ0 * (@view lpd[:, end])'
        # TODO: optimize
    end
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
    ws = make_workspace(cache, u)
    weight = cache.gauss_quadrature_weight

    phase_condition = zero(eltype(u))
    for j in 1:num_mesh(cache)
        x, dx = collocation!(ws, u, j)
        dv = reference_diff!(cache, j)
        @views for (i, t) in enumerate(turns(cache, j))
            r = diff_range(cache, j, i)
            f = H[r]
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
