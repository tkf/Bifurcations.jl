using DiffEqBase: numargs
using StaticArrays: SVector, push
using ForwardDiff


struct FixedPointBifurcationProblem{iip, HJ, H, U, T, P,
                                    } <: AbstractContinuationProblem{iip}
    homotopy_jacobian::HJ
    homotopy::H
    u0::U
    t0::T
    t_domain::Tuple{T, T}
    p::P

    # TODO: Define domain for u.  Maybe use Domains.jl?

    function FixedPointBifurcationProblem{iip}(
            homotopy::H, u0::U, t0::Real, t_domain::Tuple,
            p::P = nothing;
            homotopy_jacobian::HJ = nothing,
            ) where{iip, HJ, H, U, P}
        T = promote_type(typeof(t0), map(typeof, t_domain)...)
        new{iip, HJ, H, U, T, P}(homotopy_jacobian, homotopy,
                                 u0, t0, t_domain, p)
    end
end
const FPBPWithHJac{iip} = FixedPointBifurcationProblem{iip, <: Function}
const FPBPNoHJac{iip} = FixedPointBifurcationProblem{iip, Void}

function FixedPointBifurcationProblem(homotopy, args...; kwargs...)
    iip = numargs(homotopy) == 4
    return FixedPointBifurcationProblem{iip}(homotopy, args...; kwargs...)
end

FixedPointBifurcationProblem(homotopy, u0::Tuple, args...; kwargs...) =
    FixedPointBifurcationProblem{false}(homotopy,
                                        SVector(u0),
                                        args...; kwargs...)


struct FixedPointBifurcationCache{P, C} <: AbstractProblemCache{P}
    prob::P
    cfg::C

    FixedPointBifurcationCache(prob::P) where {P <: FPBPWithHJac} =
        new{P, Void}(prob)

    function FixedPointBifurcationCache(prob::P) where {P <: FPBPNoHJac{true}}
        x = get_u0(prob)
        y = similar(x, length(x) - 1)
        cfg = ForwardDiff.JacobianConfig((y, x) -> residual!(y, x, prob),
                                         y, x)
        return new{P, typeof(cfg)}(prob, cfg)
    end

    function FixedPointBifurcationCache(prob::P) where {P <: FPBPNoHJac{false}}
        x = get_u0(prob)
        cfg = ForwardDiff.JacobianConfig((x) -> residual!(nothing, x, prob),
                                         x)
        return new{P, typeof(cfg)}(prob, cfg)
    end
end

get_prob_cache(prob::FixedPointBifurcationProblem) =
    FixedPointBifurcationCache(prob)

function get_u0(prob::FixedPointBifurcationProblem{iip, HJ, H, U},
                ) where {iip, HJ, H, U <: AbstractVector}
    u0 = similar(prob.u0, length(prob.u0) + 1)
    u0[1:end-1] = prob.u0
    u0[end] = prob.t0
end

function get_u0(prob::FixedPointBifurcationProblem{iip, HJ, H, U},
                ) where {iip, HJ, H, U <: Real}
    return SVector(prob.u0, prob.t0)
end

function get_u0(prob::FixedPointBifurcationProblem{iip, HJ, H, U},
                ) where {iip, HJ, H, U <: SVector}
    return push(prob.u0, prob.t0)
end

residual!(H, u, cache::_C{<: FixedPointBifurcationProblem}) =
    residual!(H, u, cache.prob)

function isindomain(u, cache::_C{<: FixedPointBifurcationProblem})
    t = u[end]
    tmin, tmax = cache.prob.t_domain
    return tmin <= t <= tmax
end

# ----------------------------------------------------------- inplace interface

function residual!(H, u, prob::FixedPointBifurcationProblem{true})
    x = @view u[1:end-1]
    t = u[end]
    prob.homotopy(H, x, prob.p, t)
    return H
end

function residual_jacobian!(H, J, u, cache::_C{<: FPBPWithHJac{true}})
    prob = cache.prob
    x = @view u[1:end-1]
    t = u[end]
    prob.homotopy_jacobian(J, x, prob.p, t)
    return (H, J)
end

function residual_jacobian!(H, J, u, cache::_C{<: FPBPNoHJac{true}})
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

function residual!(_, u, prob::FixedPointBifurcationProblem{false})
    x = u[1:end-1]
    t = u[end]
    return prob.homotopy(x, prob.p, t)
end

function residual_jacobian!(_H, _J, u, cache::_C{<: FPBPWithHJac{false}})
    prob = cache.prob
    x = u[1:end-1]
    t = u[end]
    return prob.homotopy_jacobian(x, prob.p, t)
end

function residual_jacobian!(_H, _, u, cache::_C{<: FPBPNoHJac{false}})
    # TODO: No way to do this in ForwardDiff?  Or is it already
    # maximally efficient?
    H = residual!(_H, u, cache)
    J = ForwardDiff.jacobian(
        (x) -> residual!(_H, x, cache),
        u,  # x
        cache.cfg,
    )
    return (H, J)
end
