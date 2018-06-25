using ForwardDiff
using Setfield: Lens, set, get

using ...Bifurcations: maybe_subtract!

"""
Codimension-2 fixed point bifurcation problem wrapper for DifferentialEquations.
"""
struct DiffEqCodim2Problem{
        skind <: StateKind,
        tkind <: TimeKind,
        ASType <: AugmentedSystem,
        P, X, V, T, L1, L2,
        } <: Codim2Problem{skind, tkind}
    statekind::skind
    timekind::tkind
    augmented_system::ASType
    de_prob::P
    x0::X
    v0::V
    t0::T
    t_domain::Tuple{T, T}
    param_axis1::L1
    param_axis2::L2
end

StateKind(::Type{<: DiffEqCodim2Problem{skind}}) where skind = skind()
TimeKind(::Type{<: DiffEqCodim2Problem{_sk, tkind}}) where {_sk, tkind} = tkind()

function DiffEqCodim2Problem(
        de_prob,
        param_axis1::Lens,
        param_axis2::Lens,
        t_domain::Tuple;
        x0 = copy(de_prob.u0),
        v0 = copy(x0),  # TODO: fix it. this default won't work
        t0 = SVector(get(param_axis1, de_prob.p),
                     get(param_axis2, de_prob.p)),
        augmented_system = BackReferencingAS(),
        )
    t_domain = (SVector{2, eltype(t0)}(t_domain[1]),
                SVector{2, eltype(t0)}(t_domain[2]))
    return DiffEqCodim2Problem(
        statekind(de_prob), timekind(de_prob),
        augmented_system,
        de_prob, x0, v0, t0, t_domain,
        param_axis1, param_axis2)
end

AugmentedSystem(::Type{<: DiffEqCodim2Problem{_sk, _tk, ASType}}
                ) where {_sk, _tk, ASType} = ASType()

_get_augsys_cache(::BackReferencingAS, prob::DiffEqCodim2Problem) =
    BackReferencingASCache(copy(prob.v0))

# ------------------------------------------------ DiffEqCodim2BifurcationCache

struct DiffEqCodim2BifurcationCache{P, C, AC} <: AbstractProblemCache{P}
    prob::P
    augsys_cache::AC
    cfg::C
end

function DiffEqCodim2BifurcationCache(prob::DiffEqCodim2Problem)
    augsys_cache = get_augsys_cache(prob)
    return DiffEqCodim2BifurcationCache(
        prob,
        augsys_cache,
        setup_fd_config(statekind(prob), prob, augsys_cache))
end

function setup_fd_config(::MutableState, prob, augsys_cache)
    x = vcat(prob.x0, prob.v0, prob.t0)  # length = 2N + 1
    y = copy(x)
    return ForwardDiff.JacobianConfig(
        (y, x) -> _residual!(y, x, prob, augsys_cache, statekind(cache.prob)),
        y,
        x)
end

function setup_fd_config(::ImmutableState, prob, augsys_cache)
    x = vcat(prob.x0, prob.v0, prob.t0)  # length = 2N + 1
    return ForwardDiff.JacobianConfig(
        (x) -> _residual!(x, x, prob, augsys_cache, statekind(cache.prob)),
        x)
end

# ------------------------------------------------------ continuation interface

get_prob_cache(prob::DiffEqCodim2Problem) =
    DiffEqCodim2BifurcationCache(prob)

get_u0(prob::DiffEqCodim2Problem) = _get_u0(prob, prob.x0)

function _get_u0(prob::DiffEqCodim2Problem, ::AbstractVector)
    return vcat(prob.x0, as_reals(prob.v0), prob.t0)
end

function _get_u0(prob::DiffEqCodim2Problem, ::SVector) :: SVector
    return vcat(prob.x0, as_reals(prob.v0), prob.t0)
end

function _get_u0(prob::DiffEqCodim2Problem, ::Real)
    return SVector(prob.x0, as_reals(prob.v0), prob.t0...)
end

function isindomain(u, cache::DiffEqCodim2BifurcationCache)
    t = SVector(u[end-1], u[end])
    tmin, tmax = cache.prob.t_domain
    return all(tmin .<= t .<= tmax)
end

residual!(H, u, cache::DiffEqCodim2BifurcationCache) =
    _residual!(H, u, cache.prob, cache.augsys_cache, statekind(cache.prob))

residual_jacobian!(H, J, u, cache::DiffEqCodim2BifurcationCache) =
    _residual_jacobian!(H, J, u, cache, statekind(cache.prob))

# ------------------------------------------------------------------- utilities

function modified_param!(p::DiffEqCodim2Problem, u)
    t1, t2 = u[end-1:end]
    q0 = p.de_prob.p
    q1 = set(p.param_axis1, q0, t1)
    q2 = set(p.param_axis2, q1, t2)
    return q2
end

# ------------------------------------------------------------------- residual!

eigvec_constraint(v, ::NormalizingASCache) = v ⋅ v - 1
eigvec_constraint(v, augsys_cache::BackReferencingASCache) =
    augsys_cache.v ⋅ v - 1

ds_eigvec(prob::DiffEqCodim2Problem, u::AbstractArray) =
    _ds_eigvec(eltype(prob.v0), u::AbstractArray)

function _ds_eigvec(::Type{E}, u::AbstractArray) where {E <: Real}
    d = dims_from_augsys(length(u), E)
    return @view u[eigvec_range(d)]
end

function _ds_eigvec(::Type{E}, u::Array{T}) where {E <: Complex, T}
    d = dims_from_augsys(length(u), E)
    # return ComplexView(u, eigvec_range(d))  # TODO: implement
    ofs, r = divrem(d.ds_dim, 2)
    @assert r == 0
    uc = reinterpret(Complex{T}, u)
    return @view uc[ofs + 1:ofs + d.ds_dim]
end

function _ds_eigvec_expr(C::Type{<: Real}, indices)
    values = [:(u[$i]) for i in indices]
    return Expr(:call, SVector, values...)
end

function _ds_eigvec_expr(::Type{<: Complex}, indices)
    values = [:(u[$i] + u[$i + 1] * im) for i in indices[1:2:end]]
    return Expr(:call, SVector, values...)
end

# Disambiguation:
_ds_eigvec(T::Type{<: Real}, u::SVector) = _ds_eigvec_svec(T, u)
_ds_eigvec(T::Type{<: Complex}, u::SVector) = _ds_eigvec_svec(T, u)

@generated function _ds_eigvec_svec(::Type{eigvec_eltype}, u::SVector{S, T}
                                    ) where {eigvec_eltype, S, T}
    d = dims_from_augsys(S, eigvec_eltype)
    return _ds_eigvec_expr(eigvec_eltype, eigvec_range(d))
end

function _residual!(H, u, prob::DiffEqCodim2Problem,
                    augsys_cache,
                    ::MutableState)
    q = modified_param!(prob, u)

    H1, H2, H3 = output_vars(H)
    x = ds_state(u)
    v = ds_eigvec(prob, u)

    # TODO: don't allocate J
    J = ForwardDiff.jacobian(
        (dx, x) -> prob.de_prob.f(dx, x, q, 0),
        H1,  # dx
        x,
        # TODO: setup cache
    )

    A_mul_B!(H2, J, v)
    H3[] = eigvec_constraint(v, augsys_cache)

    maybe_subtract!(H1, x, statekind(prob), timekind(prob))
    maybe_subtract!(H2, v, statekind(prob), timekind(prob))
    return H
end

function _residual!(::Any, u, prob::DiffEqCodim2Problem,
                    augsys_cache,
                    ::ImmutableState)
    q = modified_param!(prob, u)

    # TODO: Can I compute H and J in one go?  Or is it already
    # maximally efficient?
    x = ds_state(u)
    H1 = prob.de_prob.f(x, q, 0)
    J = ForwardDiff.jacobian(
        (x) -> prob.de_prob.f(x, q, 0),
        x,
        # TODO: setup cache
    )

    v = ds_eigvec(prob, u)
    H2 = J * v
    H3 = eigvec_constraint(v, augsys_cache)

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
        cache.cfg,
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
        cache.cfg,
    )
    return (H, J)
end
