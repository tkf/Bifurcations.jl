module SymDiff

using Compat.Test
using Compat: Nothing
using Parameters: @with_kw, @unpack
using Setfield: set
using StaticArrays: SVector

using SymEngine: @vars, diff, subs, symbols
import SymEngine

using Bifurcations.ArrayUtils: cast_container, _normalize!
using Bifurcations.BifurcationsBase: contkind, HopfCont, SaddleNodeCont
using Bifurcations.Codim2: DiffEqCodim2Problem, get_augsys_cache, augsys,
    NormalizingAS, BackReferencingAS, ds_state, ds_eigvec
using Bifurcations.Continuations: get_prob_cache, get_u0, residual
import Bifurcations.Continuations: residual_jacobian

sjac(f, x) = [diff(f[i], x[j]) for i in 1:length(f), j in 1:length(x)]

@with_kw struct SymDiffCodim2Problem
    prob::DiffEqCodim2Problem
    param1::SymEngine.Basic
    param2::SymEngine.Basic
    x::Vector
    v::Vector
    v0::Union{Vector{SymEngine.Basic}, Nothing}
    w::Union{SymEngine.Basic, Nothing}
    f::Vector
    H::Vector
    dfdx::Matrix
    dHdu::Matrix
end

interleave(x, y) = [x y]'[:]

function SymDiffCodim2Problem(prob::DiffEqCodim2Problem)
    @vars param1 param2
    p = prob.de_prob.p
    p = set(prob.param_axis1, p, param1)
    p = set(prob.param_axis2, p, param2)

    n = length(prob.de_prob.u0)
    x = [symbols("x_$i") for i in 1:n]
    # x = cast_container(typeof(prob.de_prob.u0), x)
    x = SVector(x...)
    f = prob.de_prob.f(x, p, prob.de_prob.tspan[1])
    dfdx = sjac(f, x)

    if contkind(prob) isa SaddleNodeCont
        w = nothing
        v = [symbols("v_$i") for i in 1:n]
        if prob.augmented_system isa NormalizingAS
            v0 = v
        else
            v0 = [symbols("v0_$i") for i in 1:n]
        end
        u = vcat(x, v, param1, param2)
        H = vcat(
            f,
            dfdx * v,
            v0 ⋅ v - 1,
        )
        if prob.augmented_system isa NormalizingAS
            v0 = nothing
        end
    else
        @assert contkind(prob) isa HopfCont
        @vars w
        vr = [symbols("vr_$i") for i in 1:n]
        vi = [symbols("vi_$i") for i in 1:n]
        v = interleave(vr, vi)
        u = vcat(x, v, w, param1, param2)
        H_head = vcat(
            f,
            interleave(
                dfdx * vr .+ w .* vi,
                dfdx * vi .- w .* vr,
            ),
        )
        if prob.augmented_system isa NormalizingAS
            H = vcat(
                H_head,
                vr ⋅ vr + vi ⋅ vi - 1,
                vr ⋅ vi,
            )
            v0 = nothing
        else
            v0r = [symbols("v0r_$i") for i in 1:n]
            v0i = [symbols("v0i_$i") for i in 1:n]
            v0 = interleave(v0r, v0i)
            H = vcat(
                H_head,
                v0r ⋅ vr + v0i ⋅ vi - 1,
                v0r ⋅ vi - v0i ⋅ vr,
            )
        end
    end

    dHdu = sjac(H, u)  # aka J

    return SymDiffCodim2Problem(
        prob,
        param1,
        param2,
        x,
        v,
        v0,
        w,
        f,
        H,
        dfdx,
        dHdu,
    )
end

@with_kw struct SymDiffCodim2Cache{AC}
    prob::SymDiffCodim2Problem
    augsys_cache::AC
end

SymDiffCodim2Cache(prob::DiffEqCodim2Problem) =
    SymDiffCodim2Cache(SymDiffCodim2Problem(prob))
SymDiffCodim2Cache(prob) = SymDiffCodim2Cache(prob, get_augsys_cache(prob.prob))

function all_symbols(prob::SymDiffCodim2Problem)
    @unpack param1, param2, x, v, v0, w = prob
    syms = Any[x, v]
    v0 !== nothing && push!(syms, v0)
    w  !== nothing && push!(syms, w)
    return vcat(syms..., param1, param2)
end

all_symbols(prob_cache::SymDiffCodim2Cache) = all_symbols(prob_cache.prob)

function all_values(u, prob_cache::SymDiffCodim2Cache)
    param1, param2 = u[end-1:end]
    x = ds_state(prob_cache.prob.prob, u)
    v = ds_eigvec(prob_cache.prob.prob, u)
    if contkind(prob_cache.prob.prob) isa SaddleNodeCont
        vals = Any[x, v]
    else
        @assert contkind(prob_cache.prob.prob) isa HopfCont
        vals = Any[x, interleave(real(v), imag(v))]
    end
    if contkind(prob_cache.prob.prob) isa HopfCont
        if augsys(prob_cache.prob.prob) isa BackReferencingAS
            push!(vals, interleave(
                real(prob_cache.augsys_cache.v),
                imag(prob_cache.augsys_cache.v),
            ))
        end
        push!(vals, u[end-2])  # w
    else
        if augsys(prob_cache.prob.prob) isa BackReferencingAS
            push!(vals, prob_cache.augsys_cache.v)
        end
    end
    return vcat(vals..., param1, param2) :: AbstractVector{<: Real}
end

function evalenv(u, prob_cache::SymDiffCodim2Cache)
    syms = all_symbols(prob_cache)
    vals = all_values(u, prob_cache)
    if length(syms) != length(vals)
        error("""Number of symbols and values are different:
        #symbols = $(length(syms))
        #values  = $(length(vals))
        symbols  = $syms
        values   = $vals
        """)
    end
    return zip(syms, vals)
end

evalsym(e, env) = Float64(subs(e, env...))
evalsym(e::AbstractArray, env) = evalsym.(e, Ref(env))

function residual_jacobian(u, prob_cache::SymDiffCodim2Cache)
    env = evalenv(u, prob_cache)
    H = evalsym(prob_cache.prob.H, env)
    J = evalsym(prob_cache.prob.dHdu, env)
    return H, J
end

function _randn(rng, ::Type{Complex{T}}, n) where T
    reals = randn(rng, T, 2n)
    cs = reinterpret(Complex{T}, reals)
    @assert length(cs) == n
    return cs
end

_randn(args...) = randn(args...)

function test_residual_jacobian(prob::DiffEqCodim2Problem;
                                num = 1,
                                rng = MersenneTwister(0))
    sd_cache = SymDiffCodim2Cache(prob)
    prob_cache = get_prob_cache(prob)
    u0 = get_u0(prob)
    for _ in 1:num
        u = typeof(u0)(randn(rng, eltype(u0), length(u0)))
        if prob.augmented_system isa BackReferencingAS
            let v = sd_cache.augsys_cache.v
                v = typeof(v)(_randn(rng, eltype(v), length(v)))
                v = _normalize!(v)
                sd_cache.augsys_cache.v = v
                prob_cache.augsys_cache.v = v
            end
        end

        H1_actual = residual(u, prob_cache)
        H2_actual, J_actual = residual_jacobian(u, prob_cache)
        H_desired, J_desired = residual_jacobian(u, sd_cache)

        @test H1_actual ≈ H_desired
        @test H2_actual ≈ H_desired
        @test J_actual ≈ J_desired
    end
end

end  # module
