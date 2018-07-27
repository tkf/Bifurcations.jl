module Reparametrization

using DiffEqBase: ODEProblem, isinplace
using Parameters: @with_kw, @unpack
using Setfield: @lens, compose
using StaticArrays: SVector, SMatrix
import Setfield

using ...ArrayUtils: _zeros, badeltype
# using ...BifurcationsBase: TimeKind

struct Reparametrizer{DIM_ALL, DIM_ORIG, DIM_EXTRA, TF, TP, TM, TS, TE}
    orig_f::TF
    orig_p::TP
    M::TM
    shift::TS
    extra::TE
    invM::TM
end

function Reparametrizer(orig_f::TF, orig_p::TP, M, shift::TS, extra::TE,
                        invM = inv(M),
                        ) where {TF, TP, TS, TE}
    @assert M * invM ≈ eye(M)
    orig_TM = typeof(M)
    M, invM = promote(M, invM)
    TM = typeof(M) :: Type{<: AbstractMatrix}
    if orig_TM <: SMatrix
        @assert TM <: SMatrix
    end
    DIM_ALL = size(M, 1)
    DIM_EXTRA = length(extra)
    DIM_ORIG = DIM_ALL - DIM_EXTRA
    return Reparametrizer{DIM_ALL, DIM_ORIG, DIM_EXTRA, TF, TP, TM, TS, TE}(
        orig_f, orig_p, M, shift, extra, invM,
    )
end

function Reparametrizer{DIM_ALL, DIM_ORIG, DIM_EXTRA}(
        orig_f::TF, orig_p::TP, M::TM, shift::TS, extra::TE, invM::TM,
        ) where {DIM_ALL, DIM_ORIG, DIM_EXTRA, TF, TP, TM, TS, TE}
    return Reparametrizer{DIM_ALL, DIM_ORIG, DIM_EXTRA, TF, TP, TM, TS, TE}(
        orig_f, orig_p, M, shift, extra, invM,
    )
end

Setfield.constructor_of(::Type{<: Reparametrizer{DIM_ALL, DIM_ORIG, DIM_EXTRA}}) where {DIM_ALL, DIM_ORIG, DIM_EXTRA} =
    Reparametrizer{DIM_ALL, DIM_ORIG, DIM_EXTRA}

@generated function split_x(x::SVector{S, T},
                            rp::Reparametrizer{DIM_ALL, DIM_ORIG, DIM_EXTRA},
                            ) where {S, T, DIM_ALL, DIM_ORIG, DIM_EXTRA}
    @assert S == DIM_ALL
    idx_orig = 1:DIM_ORIG
    idx_extra = DIM_ORIG + (1:DIM_EXTRA)
    quote
        (SVector{$DIM_ORIG, $T}($((:(x[$i]) for i in idx_orig)...)),
         SVector{$DIM_EXTRA, $T}($((:(x[$i]) for i in idx_extra)...)))
    end
end

forward(x, rp::Reparametrizer) = tanh.(rp.invM * (@. x - rp.shift))
backward(y, rp::Reparametrizer) = (rp.M * atanh.(y) .+ rp.shift)

function f(y::TY, rp::Reparametrizer, t) where {TY <: SVector}
    # Original ODE:          dx/dt = f₀(x)
    # "Re-parametrization":     x = g(y) = backward(y)
    # Re-parametrized ODE:  dy/dt = dy/dx dx/dt
    #                             = (g'(y))⁻¹ f₀(g(y))
    x = backward(y, rp)
    # @assert ! badeltype(x)
    x0, x1 = split_x(x, rp)
    dx0 = rp.orig_f(x0, rp.orig_p, t)
    dx1 = rp.extra .* x1
    # @assert ! badeltype(dx1)
    dx = vcat(dx0, dx1) :: SVector

    # (g'(y))⁻¹ = (g⁻¹(x))' = (tanh(M⁻¹ (x - shift)))' = tanh' M⁻¹
    dy = (1 .- y .* y) .* (rp.invM * dx)
    dy = SVector(dy...)  # TODO: don't
    @assert ! badeltype(dy)
    return dy
    # return dy :: TY
end
# Cannot use `dy :: TY` because of the type instability in
# first_lyapunov_coefficient.

f(u::A, p::Reparametrizer{DIM_ALL}, t) where {A <: AbstractVector, DIM_ALL} =
    A(f(SVector{DIM_ALL}(u), p, t))

function f(du, u, rp::Reparametrizer{DIM_ALL}, t) where {DIM_ALL}
    du .= f(SVector{DIM_ALL}(u), rp, t)
    nothing
end

function reparametrize(
        ode::ODEProblem, M::AbstractMatrix;
        shift = 0.0,
        extra::Union{AbstractVector, Real} = SVector{0, Float64}(),
        )
    rp = Reparametrizer(ode.f, ode.p, M, shift, extra)
    x1 = _zeros(ode.u0, length(extra))
    x = vcat(ode.u0, x1)
    y = forward(x, rp)
    u0 = typeof(x)(y)  # since `forward` returns `SVector`
    return ODEProblem{isinplace(ode)}(f, u0, ode.tspan, rp)
end

function reparametrize(
        ode::ODEProblem, rng::AbstractRNG;
        orthogonal = true,
        extra::Union{AbstractVector, Real} = SVector{0, Float64}(),
        kwargs...)
    tries = 10
    dim = length(ode.u0) + length(extra)
    for _ in 1:tries
        M = SMatrix{dim, dim}(randn(rng, dim, dim))
        if orthogonal
            M, _ = qr(M)
        end
        try
            invM = inv(M)
        catch
            continue
        end
        return reparametrize(ode, M; extra=extra, kwargs...)
    end
    error("Cannot obtain an invertible matrix after $tries tries.")
end

reparametrize(ode::ODEProblem; seed::Integer = 0, kwargs...) =
    reparametrize(ode, MersenneTwister(seed); kwargs...)


using ...Bifurcations: FixedPointBifurcationProblem, BifurcationProblem,
    unbounded_phase_space

const orig_p = @lens _.orig_p

function reparametrize(prob::FixedPointBifurcationProblem,
                       args...; kwargs...)

    let u0 = prob.p.de_prob.u0
        # TODO: support transforming prob.u0 and prob.phase_space:
        @assert prob.u0 ≈ u0
        @assert prob.phase_space == unbounded_phase_space(u0)
    end
    @assert prob.homotopy_jacobian == nothing

    ode = reparametrize(prob.p.de_prob, args...; kwargs...)
    return BifurcationProblem(
        ode,
        compose(prob.p.param_axis, orig_p),
        prob.t_domain,
    )
end

end  # module
