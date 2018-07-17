using ForwardDiff
using Parameters: @with_kw, @unpack
using Setfield: set
using StaticArrays: SVector

using ..CompatUtils: @required
import ..Continuations: get_prob_cache, get_u0, residual!, residual_jacobian!,
    isindomain

using ..Codim1LimitCycle
using ..Codim1LimitCycle: LimitCycleProblem, LimitCycleCache, residual_lc!

abstract type
    Codim2LCProblem{skind, tkind} <: BifurcationProblem{skind, tkind}
end


struct FoldLimitCycleProblem{
        skind <: MutableState,  # <: StateKind,
        tkind <: TimeKind,
        P <: LimitCycleProblem,
        XS <: AbstractMatrix,
        R <: Real,
        T, L1, L2,
        } <: Codim2LCProblem{skind, tkind}
    statekind::skind
    timekind::tkind
    super::P
    xs0::XS
    vs0::XS
    l0::R
    dl0::R
    t0::T
    t_domain::Tuple{T, T}
    param_axis1::L1
    param_axis2::L2
end

function FoldLimitCycleProblem(
        prob::LimitCycleProblem;
        (@required vs0),
        (@required dl0),
        xs0 = prob.xs0,
        l0 = prob.l0,
        param_axis1 = prob.param_axis1,
        param_axis2 = prob.param_axis2,
        t0 = SVector(get(param_axis1, prob.de_prob.p),
                     get(param_axis2, prob.de_prob.p)),
        t2_domain = (typemin(eltype(t0)), typemax(eltype(t0))),
        t_domain = (
            SVector(prob.t_domain[1], t2_domain[1]),
            SVector(prob.t_domain[2], t2_domain[2]),
        ),
        )
    return FoldLimitCycleProblem(
        MutableState(),
        timekind(prob.de_prob),
        prob,
        xs0,
        vs0,
        l0,
        dl0,
        t0,
        t_domain,
        param_axis1,
        param_axis2,
    )
end

function FoldLimitCycleProblem(;
        (@required vs0),
        (@required t0),
        (@required param_axis1),
        (@required param_axis2),
        kwargs...)
    prob = LimitCycleProblem(;
        t0 = t0[1],
        param_axis = param_axis1,
        kwargs...)
    return FoldLimitCycleProblem(
        prob;
        vs0 = vs0,
        t0 = t0,
        param_axis1 = param_axis1,
        param_axis2 = param_axis2,
    )
end

# ----------------------------------------------------------------------- cache

struct FoldLimitCycleCache{
        P <: FoldLimitCycleProblem,
        C <: LimitCycleCache,
        } <: AbstractProblemCache{P}
    super::C
    prob::P
end

function FoldLimitCycleCache(prob)
    return FoldLimitCycleCache(
        LimitCycleCache(prob.super),
        prob,
    )
end

# ------------------------------------------------------ continuation interface

get_prob_cache(prob::FoldLimitCycleProblem) = FoldLimitCycleCache(prob)

get_u0(prob::FoldLimitCycleProblem) = vcat(
    view(prob.xs0, :),          # \__ ξ
    prob.l0,                    # /
    view(prob.vs0, :),          # \__ η
    prob.dl0,                   # /
    prob.t0,                    # --- α
)
# modified_param! requires prob.t0 to be at the end

Codim1LimitCycle.u_idx_param(cache::FoldLimitCycleCache) =
    2(length(cache.prob.xs0) + 1) + (1:2)

isindomain(u, cache::FoldLimitCycleCache) =
    Codim1LimitCycle.isindomain_lc(u, cache)

# ------------------------------------------------------------------- residual!

# H = [
#     H_lc = [                                           # F(ξ, α) = 0
#         state_condition
#         phase_condition
#     ]
#     eigvec_condition                                   # ∂₁ F(ξ, α) η = 0
#     eigvec_constraint                                  # <η,η> - 1 = 0
# ]

# Almost identical to: [[../codim2/diffeq.jl::modified_param!]]
function modified_param!(p, u)
    t1, t2 = u[end-1:end]
    q0 = p.super.de_prob.p
    q1 = set(p.param_axis1, q0, t1)
    q2 = set(p.param_axis2, q1, t2)
    return q2
end

function residual!(H, u, cache::FoldLimitCycleCache)
    prob = cache.prob
    q = modified_param!(prob, u)

    n_lc = length(prob.xs0) + 1
    ξ = view(u, 1:n_lc)                  # [xs; period]
    η = view(u, (1:n_lc) + n_lc)         # ∂ξ
    H_lc = view(H, 1:n_lc)               # F(ξ, α)
    H_ev = view(H, (1:n_lc) + n_lc)      # ∂₁ F(ξ, α) η
    H_ec = @view H[end]                  # <η,η> - 1

    ForwardDiff.jacobian!(
        reshape(H_ev, :, 1),  # =  ∂₁ F(ξ, α) η  =  ∂ᵧ F(ξ + γ η, α)
        (H_lc, d) -> residual_lc!(H_lc, (@. ξ + d[1] * η), q, cache.super),
        H_lc,
        SVector(zero(eltype(η))),
        # TODO: setup cache
    )

    H_ec[] = η ⋅ η - 1

    return H
end

# ---------------------------------------------------------- residual_jacobian!

function residual_jacobian!(H, J, u, cache::FoldLimitCycleCache)
    # TODO: implement residual_jacobian! manually
    ForwardDiff.jacobian!(
        J,
        (y, x) -> residual!(y, x, cache),
        H,  # y
        u,  # x
    )
    return (H, J)
end
