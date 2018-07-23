using ForwardDiff
using Setfield: Lens, set, get
using StaticArrays: @SMatrix, SMatrix, SVector

using ...Bifurcations: maybe_subtract!

"""
Codimension-2 fixed point bifurcation problem wrapper for DifferentialEquations.
"""
struct DiffEqCodim2Problem{
        skind <: StateKind,
        tkind <: TimeKind,
        ASType <: AugmentedSystem,
        P, X, V, W, T, L1, L2,
        } <: Codim2Problem{skind, tkind}
    statekind::skind
    timekind::tkind
    augmented_system::ASType
    de_prob::P
    x0::X
    v0::V
    w0::W
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
        w0 = 0,
        t0 = SVector(get(param_axis1, de_prob.p),
                     get(param_axis2, de_prob.p)),
        augmented_system = BackReferencingAS(),
        )
    t_domain = (SVector{2, eltype(t0)}(t_domain[1]),
                SVector{2, eltype(t0)}(t_domain[2]))
    return DiffEqCodim2Problem(
        statekind(de_prob), timekind(de_prob),
        augmented_system,
        de_prob, x0, v0, w0, t0, t_domain,
        param_axis1, param_axis2)
end

AugmentedSystem(::Type{<: DiffEqCodim2Problem{_sk, _tk, ASType}}
                ) where {_sk, _tk, ASType} = ASType()

_get_augsys_cache(::BackReferencingAS, prob::DiffEqCodim2Problem) =
    BackReferencingASCache(copy(prob.v0))

# ------------------------------------------------ DiffEqCodim2BifurcationCache

struct DiffEqCodim2BifurcationCache{P, C, AC, F} <: AbstractProblemCache{P}
    prob::P
    augsys_cache::AC
    cfg::C
    residual::F
end

function DiffEqCodim2BifurcationCache(prob::DiffEqCodim2Problem)
    augsys_cache = get_augsys_cache(prob)
    return _DiffEqCodim2BifurcationCache(statekind(prob), prob, augsys_cache)
end

function _DiffEqCodim2BifurcationCache(::MutableState, prob, augsys_cache)
    x = copy(get_u0(prob))
    y = similar(x, length(x) - 1)
    residual = (y, x) -> _residual!(y, x, prob, augsys_cache)
    cfg = ForwardDiff.JacobianConfig(residual, y, x)
    return DiffEqCodim2BifurcationCache(
        prob,
        augsys_cache,
        cfg,
        residual)
end

function _DiffEqCodim2BifurcationCache(::ImmutableState, prob, augsys_cache)
    x = copy(get_u0(prob))
    residual = (x) -> _residual!(nothing, x, prob, augsys_cache)
    cfg = ForwardDiff.JacobianConfig(residual, x)
    return DiffEqCodim2BifurcationCache(
        prob,
        augsys_cache,
        cfg,
        residual)
end

# ------------------------------------------------------ continuation interface

get_prob_cache(prob::DiffEqCodim2Problem) =
    DiffEqCodim2BifurcationCache(prob)

get_u0(prob::DiffEqCodim2Problem) = _get_u0(prob, prob.x0)

function _get_u0(prob::DiffEqCodim2Problem, ::AbstractVector)
    return vcat(prob.x0, as_reals(prob.v0), imag_eigval(prob), prob.t0)
end

function _get_u0(prob::DiffEqCodim2Problem, ::SVector) :: SVector
    return vcat(prob.x0, as_reals(prob.v0), imag_eigval(prob), prob.t0)
end

function _get_u0(prob::DiffEqCodim2Problem, ::Real)
    return SVector(prob.x0, as_reals(prob.v0)...,
                   imag_eigval(prob)..., prob.t0...)
end

imag_eigval(prob) = _imag_eigval(eltype(prob.v0), prob.w0)
_imag_eigval(::Type{T}, ::Real) where {T <: Real} = SVector{0, T}()
_imag_eigval(::Type{<: Complex}, w0::Real) = SVector(w0)

isindomain(u, cache::DiffEqCodim2BifurcationCache) =
    _isindomain(u, cache, eltype(cache.prob.v0))

function _isindomain(u, cache, ::Type{<: Real})
    t = SVector(u[end-1], u[end])
    tmin, tmax = cache.prob.t_domain
    return all(tmin .<= t .<= tmax)
end

function _isindomain(u, cache, ::Type{<: Complex})
    # Constrain to positive angular velocity (w = u[end-2])
    p = SVector(u[end-2], u[end-1], u[end])
    tmin, tmax = cache.prob.t_domain
    pmin = SVector(zero(tmin[1]), tmin[1], tmin[2])
    pmax = SVector(typemax(tmax[1]), tmax[1], tmax[2])
    return all(pmin .<= p .<= pmax)
end

residual!(H, u, cache::DiffEqCodim2BifurcationCache) =
    _residual!(H, u, cache.prob, cache.augsys_cache)

_residual!(H, u, prob, augsys_cache) =
    _residual!(H, u, prob, augsys_cache, statekind(prob))

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

eigvec_constraint(v::AbstractVector{T}, ::NormalizingASCache) where {T <: Real} =
    SVector{1, T}(v ⋅ v - 1)

eigvec_constraint(v::AbstractVector{<: Complex{T}}, ::NormalizingASCache) where {T} =
    SVector{2, T}(real(v ⋅ v) - 1, real(v) ⋅ imag(v))

eigvec_constraint(v::AbstractMatrix{T}, ::NormalizingASCache) where {T} =
    SVector{2, T}(
        v[:] ⋅ v[:] - 1,
        v[:, 1] ⋅ v[:, 2],
    )

eigvec_constraint(v::AbstractVector, augsys_cache::BackReferencingASCache) =
    as_reals(augsys_cache.v ⋅ v - 1)

function eigvec_constraint(v::AbstractMatrix{T},
                           augsys_cache::BackReferencingASCache,
                           ) where {T}
    vr = @view v[:, 1]
    vi = @view v[:, 2]
    v0r = real(augsys_cache.v)
    v0i = imag(augsys_cache.v)
    return SVector{2, T}(
        v0r ⋅ vr + v0i ⋅ vi - 1,
        v0r ⋅ vi - v0i ⋅ vr,
    )
end

ds_state(prob::DiffEqCodim2Problem, u::AbstractArray) =
    _ds_state(eltype(prob.v0), u)

function _ds_state(E::Type, u)
    d = dims_from_augsys(length(u), E)
    return @view u[1:d.ds_dim]
end

@generated function _ds_state(::Type{E}, u::SVector{S, T}) where {E, S, T}
    d = dims_from_augsys(S, E)
    values = [:(u[$i]) for i in 1:d.ds_dim]
    quote
        SVector{$(length(values)), $T}($(values...))
    end
end

ds_eigvecmat(prob::DiffEqCodim2Problem, u::AbstractArray) =
    _ds_eigvecmat(eltype(prob.v0), u)

_ds_eigvecmat(E::Type{<: Real}, u) = _ds_eigvec(E, u)
function _ds_eigvecmat(E::Type{<: Complex}, u)
    d = dims_from_augsys(length(u), E)
    return reshape((@view u[eigvec_range(d)]), (2, :))'  # TODO: optimize
end

@generated function _ds_eigvecmat(::Type{E}, u::SVector{S}) where {E <: Complex, S}
    d = dims_from_augsys(S, E)
    indices = eigvec_range(d)
    values = [:(u[$i]) for i in indices]
    quote
        SMatrix{2, $(length(values) ÷ 2)}($(values...))'
    end
end

ds_eigvec(prob::DiffEqCodim2Problem, u::AbstractArray) =
    _ds_eigvec(eltype(prob.v0), u::AbstractArray)

ds_eigvec(::HopfCont, u::AbstractArray) = _ds_eigvec(Complex, u)
ds_eigvec(::SaddleNodeCont, u::AbstractArray) = _ds_eigvec(Real, u)

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

ds_eigval(prob::DiffEqCodim2Problem, u::AbstractArray) =
    _ds_eigval(eltype(prob.v0), u::AbstractArray)

_ds_eigval(::Type{<: Real}, ::AbstractArray) = 0

function _ds_eigval(T::Type{<: Complex}, u::AbstractArray)
    d = dims_from_augsys(length(u), T)
    return u[eigval_index(d)] * im
end

ds_eigvalmat(prob::DiffEqCodim2Problem, u::AbstractArray) =
    _ds_eigvalmat(eltype(prob.v0), u::AbstractArray)

_ds_eigvalmat(::Type{<: Real}, ::AbstractArray) = 0

function _ds_eigvalmat(E::Type{<: Complex}, u::AbstractArray{T}) where {T}
    d = dims_from_augsys(length(u), E)
    w = u[eigval_index(d)]
    return @SMatrix T[
         0 w
        -w 0
    ]
end

state_cond_view(d::VarDims, H) = @view H[1:d.ds_dim]       # f(x)
eigvec_cond_view(d::VarDims, H) = @view H[eigvec_range(d)] # Jv - 1
eigvec_cons_view(d::VarDims, H) = @view H[last(eigvec_range(d)) + 1:end]

J_mul_v!(Jv, fx, x, v::AbstractVector{<: Real}, q, f) =
    ForwardDiff.jacobian!(
        reshape(view(Jv, :), :, 1),
        (fx, d) -> f(fx, (@. x + d[1] * v), q, 0),
        fx,
        SVector(zero(eltype(x))),
        # TODO: setup cache
    )

function J_mul_v!(Jv, fx, x, v::AbstractVector{<: Complex}, q, f)
    vr = real(v)
    vi = imag(v)
    JvT = reshape(view(Jv, :), 2, :)'  # not a view in Julia 0.6
    ForwardDiff.jacobian!(
        JvT,
        (fx, d) -> f(fx, (@. x + d[1] * vr + d[2] * vi), q, 0),
        fx,
        SVector(zero(eltype(x)), zero(eltype(x))),
        # TODO: setup cache
    )
    Jv .= view(JvT', :)
end
# It's a bit tricky to handle complex vector case (see also
# `H2[2:2:end]` below).  Maybe `reinterpret` -compatible memory layout
# was not the best choice?

function _residual!(H, u, prob::DiffEqCodim2Problem,
                    augsys_cache,
                    ::MutableState)
    q = modified_param!(prob, u)

    dim = dims_from_augsys(length(u), contkind(prob))
    H1 = state_cond_view(dim, H)
    H2 = eigvec_cond_view(dim, H)
    H3 = eigvec_cons_view(dim, H)
    x = ds_state(prob, u)
    v = ds_eigvec(prob, u)
    iw = ds_eigval(prob, u)

    # TODO: don't allocate x + d * v
    J_mul_v!(H2, H1, x, v, q, prob.de_prob.f)
    if iw != 0
        # "H2" -= iw * v
        w = imag(iw)
        @. H2[1:2:end] += w * imag(v)
        @. H2[2:2:end] -= w * real(v)
    end
    H3[:] = eigvec_constraint(v, augsys_cache)

    maybe_subtract!(H1, x, statekind(prob), timekind(prob))
    maybe_subtract!(H2, v, statekind(prob), timekind(prob))
    return H
end

J_mul_v(x, v::AbstractVector{<: Real}, q, f) =
    ForwardDiff.jacobian(
        (d) -> f((@. x + d[1] * v), q, 0),
        SVector(zero(eltype(x))),
        # TODO: setup cache
    )[:]

function J_mul_v(x, v::AbstractVector{<: Complex}, q, f)
    vr = real(v)
    vi = imag(v)
    Jv = ForwardDiff.jacobian(
        (d) -> f((@. x + d[1] * vr + d[2] * vi), q, 0),
        SVector(zero(eltype(x)), zero(eltype(x))),
        # TODO: setup cache
    )
    return @. Jv[:, 1] + im * Jv[:, 2]
end

function J_mul_v(x, v::AbstractMatrix, q, f)
    vr = v[:, 1]
    vi = v[:, 2]
    Jv = ForwardDiff.jacobian(
        (d) -> f((@. x + d[1] * vr + d[2] * vi), q, 0),
        SVector(zero(eltype(x)), zero(eltype(x))),
        # TODO: setup cache
    )
    return Jv
end

function _residual!(::Any, u, prob::DiffEqCodim2Problem,
                    augsys_cache,
                    ::ImmutableState)
    q = modified_param!(prob, u)

    x = ds_state(prob, u)
    v = ds_eigvecmat(prob, u)
    iw = ds_eigvalmat(prob, u)

    # TODO: Can I compute H and Jv in one go?  Or is it already
    # maximally efficient?
    H1 = prob.de_prob.f(x, q, 0)
    Jv = J_mul_v(x, v, q, prob.de_prob.f)
    H2 = as_reals(Jv .- v * iw)

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
        cache.residual,
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
    H = cache.residual(u)
    J = ForwardDiff.jacobian(
        cache.residual,
        u,  # x
        cache.cfg,
    )
    return (H, J)
end
