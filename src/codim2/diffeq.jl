using ForwardDiff
using Setfield: Lens, set, get

using ...Bifurcations: maybe_subtract!
using ..Continuations: _similar

"""
Codimension-2 fixed point bifurcation problem wrapper for DifferentialEquations.
"""
struct DiffEqCodim2BifurcationProblem{
        skind <: StateKind,
        tkind <: TimeKind,
        P, X, V, T, L1, L2,
        } <: Codim2BifurcationProblem{skind, tkind}
    statekind::skind
    timekind::tkind
    de_prob::P
    x0::X
    v0::V
    t0::T
    t_domain::Tuple{T, T}
    param_axis1::L1
    param_axis2::L2
end

StateKind(::Type{<: DiffEqCodim2BifurcationProblem{skind}}) where skind = skind()
TimeKind(::Type{<: DiffEqCodim2BifurcationProblem{_sk, tkind}}) where {_sk, tkind} = tkind()

function DiffEqCodim2BifurcationProblem(
        de_prob,
        param_axis1::Lens,
        param_axis2::Lens,
        t_domain::Tuple;
        x0 = copy(de_prob.u0),
        v0 = copy(x0),  # TODO: fix it. this default won't work
        t0 = SVector(get(param_axis1, de_prob.p),
                     get(param_axis2, de_prob.p)),
        )
    t_domain = (SVector{2, eltype(t0)}(t_domain[1]),
                SVector{2, eltype(t0)}(t_domain[2]))
    return DiffEqCodim2BifurcationProblem(
        statekind(de_prob), timekind(de_prob),
        de_prob, x0, v0, t0, t_domain,
        param_axis1, param_axis2)
end

# ------------------------------------------------ DiffEqCodim2BifurcationCache

struct DiffEqCodim2BifurcationCache{P, C} <: AbstractProblemCache{P}
    prob::P
    cfg::C
end

DiffEqCodim2BifurcationCache(prob::DiffEqCodim2BifurcationProblem) =
    DiffEqCodim2BifurcationCache(
        prob,
        setup_fd_config(statekind(prob), prob.de_prob))

function setup_fd_config(::MutableState, de_prob)
    u0 = de_prob.u0
    y = similar(u0, length(u0) * 2 + 1)
    x = similar(u0, length(u0) * 2 + 2)
    return ForwardDiff.JacobianConfig(
        (y, x) -> de_prob.f(y, x, de_prob, 0),
        y,
        x)
end

function setup_fd_config(::ImmutableState, de_prob)
    u0 = de_prob.u0
    x = _similar(u0, length(u0) * 2 + 2)
    return ForwardDiff.JacobianConfig(
        (x) -> de_prob.f(x, de_prob, 0),
        x)
end

# ------------------------------------------------------ continuation interface

get_prob_cache(prob::DiffEqCodim2BifurcationProblem) =
    DiffEqCodim2BifurcationCache(prob)

get_u0(prob::DiffEqCodim2BifurcationProblem) = _get_u0(prob, prob.x0)

function _get_u0(prob::DiffEqCodim2BifurcationProblem, ::AbstractVector)
    return vcat(prob.x0, prob.v0, prob.t0)
end

function _get_u0(prob::DiffEqCodim2BifurcationProblem, ::SVector) :: SVector
    return vcat(prob.x0, prob.v0, prob.t0)
end

function _get_u0(prob::DiffEqCodim2BifurcationProblem, ::Real)
    return SVector(prob.x0, prob.v0, prob.t0...)
end

function isindomain(u, cache::DiffEqCodim2BifurcationCache)
    t = SVector(u[end-1], u[end])
    tmin, tmax = cache.prob.t_domain
    return all(tmin .<= t .<= tmax)
end

residual!(H, u, cache::DiffEqCodim2BifurcationCache) =
    _residual!(H, u, cache, statekind(cache.prob))

residual_jacobian!(H, J, u, cache::DiffEqCodim2BifurcationCache) =
    _residual_jacobian!(H, J, u, cache, statekind(cache.prob))

# ------------------------------------------------------------------- utilities

function modified_param!(p::DiffEqCodim2BifurcationProblem, u)
    t1, t2 = u[end-1:end]
    q0 = p.de_prob.p
    q1 = set(p.param_axis1, q0, t1)
    q2 = set(p.param_axis2, q1, t2)
    return q2
end

# ------------------------------------------------------------------- residual!

function _residual!(H, u, cache::DiffEqCodim2BifurcationCache,
                    ::MutableState)
    prob = cache.prob
    q = modified_param!(prob, u)

    H1, H2, H3 = output_vars(H)
    x = ds_state(u)
    v = ds_eigvec(u)

    # TODO: don't allocate J
    J = ForwardDiff.jacobian(
        (dx, x) -> prob.de_prob.f(dx, x, q, 0),
        H1,  # dx
        x,
        # cache.cfg,  # TODO: make it work
    )

    A_mul_B!(H2, J, v)
    H3[] = v ⋅ v - 1

    maybe_subtract!(H1, x, statekind(prob), timekind(prob))
    maybe_subtract!(H2, v, statekind(prob), timekind(prob))
    return H
end

function _residual!(::Any, u, cache::DiffEqCodim2BifurcationCache,
                    ::ImmutableState)
    prob = cache.prob
    q = modified_param!(prob, u)

    # TODO: Can I compute H and J in one go?  Or is it already
    # maximally efficient?
    x = ds_state(u)
    H1 = prob.de_prob.f(x, q, 0)
    J = ForwardDiff.jacobian(
        (x) -> prob.de_prob.f(x, q, 0),
        x,
        # cache.cfg,  # TODO: make it work
    )

    v = ds_eigvec(u)
    H2 = J * v
    H3 = v ⋅ v - 1

    return cat_outputs(
        maybe_subtract!(H1, x, statekind(prob), timekind(prob)),
        maybe_subtract!(H2, v, statekind(prob), timekind(prob)),
        H3)
end

# ---------------------------------------------------------- residual_jacobian!

function _residual_jacobian!(H, J, u, cache::DiffEqCodim2BifurcationCache,
                             ::MutableState)
    ForwardDiff.jacobian!(
        J,
        (y, x) -> residual!(y, x, cache),
        H,  # y
        u,  # x
        # TODO: setup cache
    )
    return (H, J)
end

function _residual_jacobian!(_H, _J, u, cache::DiffEqCodim2BifurcationCache,
                             ::ImmutableState)
    # TODO: Can I compute H and J in one go?  Or is it already
    # maximally efficient?
    H = residual!(_H, u, cache)
    J = ForwardDiff.jacobian(
        (x) -> residual!(_H, x, cache),
        u,  # x
        # TODO: setup cache
    )
    return (H, J)
end
