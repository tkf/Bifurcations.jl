using DiffEqBase: numargs
using StaticArrays: SVector, push
using ForwardDiff

using .Codim1: Codim1Problem


"""
Fixed point bifurcation problem.

See also: [`AbstractContinuationProblem`](@ref)

# Fields
* `homotopy::Function`:
  A function to compute ``H(x, t)`` where ``H`` is a homotopy
  ``H: \\mathbb R^N \\times \\mathbb R \\to \\mathbb R^N``.
  Function `homotopy` must be callable in one of the following form:
  `homotopy(x, p, t) ↦ H` (return `H`) for mutable state type or
  `homotopy(H, x, p, t)` (mutate `H`) for immutable state type.
* `homotopy_jacobian::Union{Function, Nothing}`:
  A function to compute ``H(x, t)`` and its Jacobian
  ``J = \\partial H / \\partial (x, t) \\in \\mathbb R^{N \\times (N+1)}``.
  Function `homotopy_jacobian` must be callable in one of the following form:
  `homotopy_jacobian(x, p, t) ↦ (H, J)` (return `(H, J)`) or
  `homotopy_jacobian(H, J, x, p, t)` (mutate `H` and `J`).
* `u0::Union{AbstractArray, Real}`: Initial state.
* `t0::Real`: Initial parameter.
* `t_domain::Tuple{<:Real, <:Real}`: Range of the parameter.
* `phase_space::Tuple{typeof(u0), typeof(u0)}`: A pair of lower and
  upper bound of the phase space.  Default is unbounded.
* `p`: Model parameter (constants).
"""
struct FixedPointBifurcationProblem{skind <: StateKind,
                                    tkind <: TimeKind,
                                    HJ, H, U, T, P,
                                    } <: Codim1Problem{skind, tkind}
    statekind::skind
    timekind::tkind
    homotopy_jacobian::HJ
    homotopy::H
    # TODO: rename u0 -> x0 (consistently use x for dynamical system state)
    u0::U
    t0::T
    t_domain::Tuple{T, T}
    phase_space::Tuple{U, U}
    p::P

    # TODO: Use Domains.jl?
end

"""
    unbounded_phase_space(u0::AbstractVector)
    unbounded_phase_space(u0::Real)
"""
function unbounded_phase_space(x0::X) where {X <: AbstractVector}
    return (X(fill(-Inf, length(x0))),
            X(fill(+Inf, length(x0))))
end
unbounded_phase_space(::X) where {X <: Real} = (X(-Inf), X(Inf))

function FixedPointBifurcationProblem(
        statekind::StateKind, timekind::TimeKind,
        homotopy::H, u0::U, t0::Real, t_domain::Tuple,
        p::P = nothing;
        homotopy_jacobian::HJ = nothing,
        phase_space = unbounded_phase_space(u0),
        ) where{HJ, H, U, P}
    T = promote_type(typeof(t0), map(typeof, t_domain)...)
    return FixedPointBifurcationProblem(
        statekind, timekind,
        homotopy_jacobian, homotopy,
        u0, t0, t_domain, phase_space, p)
end

struct HasJac end
struct NoJac end
hasjac(::FixedPointBifurcationProblem{<:Any, <:Any, <:Function}) = HasJac()
hasjac(::FixedPointBifurcationProblem{<:Any, <:Any, Void}) = NoJac()

function FixedPointBifurcationProblem(tkind::TimeKind,
                                      homotopy, args...; kwargs...)
    skind = numargs(homotopy) == 4 ? MutableState() : ImmutableState()
    return FixedPointBifurcationProblem(
        skind, tkind,
        homotopy, args...; kwargs...)
end

as_immutable_state(x::Tuple) = SVector(x)
as_immutable_state(x::Number) = x

function FixedPointBifurcationProblem(tkind::TimeKind,
                                      homotopy,
                                      u0::Union{Tuple, Number},
                                      args...; kwargs...)
    return FixedPointBifurcationProblem{false, tkind}(
        ImmutableState(), tkind,
        homotopy,
        as_immutable_state(u0),
        args...; kwargs...)
end


# TODO: trait
BifurcationsBase.contkind(::FixedPointBifurcationProblem) =
    BifurcationsBase.FixedPointCont()


struct FixedPointBifurcationCache{P, C} <: AbstractProblemCache{P}
    prob::P
    cfg::C
end

FixedPointBifurcationCache(prob::FixedPointBifurcationProblem) =
    _FixedPointBifurcationCache(statekind(prob), hasjac(prob), prob)

_FixedPointBifurcationCache(::Any, ::HasJac, prob) =
    FixedPointBifurcationCache(prob, nothing)

function _FixedPointBifurcationCache(::MutableState, ::NoJac, prob)
    x = get_u0(prob)
    y = similar(x, length(x) - 1)
    cfg = ForwardDiff.JacobianConfig((y, x) -> residual!(y, x, prob),
                                     y, x)
    return FixedPointBifurcationCache(prob, cfg)
end

function _FixedPointBifurcationCache(::ImmutableState, ::NoJac, prob)
    x = get_u0(prob)
    cfg = ForwardDiff.JacobianConfig((x) -> residual!(nothing, x, prob),
                                     x)
    return FixedPointBifurcationCache(prob, cfg)
end

TimeKind(::Type{<: FixedPointBifurcationCache{P}}) where P = TimeKind(P)

get_prob_cache(prob::FixedPointBifurcationProblem) =
    FixedPointBifurcationCache(prob)

get_u0(prob::FixedPointBifurcationProblem) = _get_u0(prob, prob.u0)

function _get_u0(prob::FixedPointBifurcationProblem, ::AbstractVector)
    u0 = similar(prob.u0, length(prob.u0) + 1)
    u0[1:end-1] = prob.u0
    u0[end] = prob.t0
    return u0
end

function _get_u0(prob::FixedPointBifurcationProblem, ::Real)
    return SVector(prob.u0, prob.t0)
end

function _get_u0(prob::FixedPointBifurcationProblem, ::SVector)
    return push(prob.u0, prob.t0)
end

residual!(H, u, cache::FixedPointBifurcationCache) =
    _residual!(H, u, cache.prob,
               statekind(cache.prob),
               cache.prob.u0)

residual_jacobian!(H, J, u, cache::FixedPointBifurcationCache) =
    _residual_jacobian!(H, J, u, cache,
                        statekind(cache.prob),
                        hasjac(cache.prob))

function isindomain(u, cache::FixedPointBifurcationCache)
    tmin, tmax = cache.prob.t_domain
    xmin, xmax = cache.prob.phase_space
    umin = vcat(xmin, SVector(tmin))
    umax = vcat(xmax, SVector(tmax))
    return all(umin .<= u .<= umax)
end

# ----------------------------------------------------------- inplace interface

function _residual!(H, u, prob::FixedPointBifurcationProblem,
                    ::MutableState, ::Any)
    x = @view u[1:end-1]
    t = u[end]
    prob.homotopy(H, x, prob.p, t)
    return H
end

function _residual_jacobian!(H, J, u, cache::FixedPointBifurcationCache,
                             ::MutableState, ::HasJac)
    prob = cache.prob
    x = @view u[1:end-1]
    t = u[end]
    prob.homotopy_jacobian(H, J, x, prob.p, t)
    return (H, J)
end

function _residual_jacobian!(H, J, u, cache::FixedPointBifurcationCache,
                             ::MutableState, ::NoJac)
    ForwardDiff.jacobian!(
        J,
        (y, x) -> residual!(y, x, cache),
        H,  # y
        u,  # x
        cache.cfg,
    )
    return (H, J)
end

# ------------------------------------------------------ out-of-place interface

function _residual!(::Any, u, prob::FixedPointBifurcationProblem,
                    ::ImmutableState, ::Real)
    x = u[1]
    t = u[end]
    return SVector(prob.homotopy(x, prob.p, t))
end

function _residual!(::Any, u, prob::FixedPointBifurcationProblem,
                    ::ImmutableState, ::AbstractArray)
    x = u[1:end-1]
    t = u[end]
    return prob.homotopy(x, prob.p, t)
end

function _residual_jacobian!(_H, _J, u, cache::FixedPointBifurcationCache,
                             ::ImmutableState, ::HasJac)
    prob = cache.prob
    x = u[1:end-1]
    t = u[end]
    return prob.homotopy_jacobian(x, prob.p, t)
end

function _residual_jacobian!(_H, _, u, cache::FixedPointBifurcationCache,
                             ::ImmutableState, ::NoJac)
    # TODO: Can I compute H and J in one go?  Or is it already
    # maximally efficient?
    H = residual!(_H, u, cache)
    J = ForwardDiff.jacobian(
        (x) -> residual!(_H, x, cache),
        u,  # x
        cache.cfg,
    )
    return (H, J)
end
