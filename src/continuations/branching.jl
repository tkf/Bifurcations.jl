using ForwardDiff
using StaticArrays: SMatrix, SArray
using Setfield: @set


"""
Find vectors tJ2 and cotJ such that ker(J) = span{tJ1, tJ2} and
ker(J') = span{cotJ} (left nullspace or cokernel):
"""
function find_more_nullspaces(Q, L, rtol, atol, max_steps)
    y = _zeros(Q, size(Q, 1) - 1)
    y = @set y[length(y)] = 1   # y[end] doesn't work

    if abs(L[end-1, end-1]) < atol
        tJ2 = Q[:, end-1]
        cotJ = _A_ldiv_B!((@view Q[1:end-1, 1:end-1]), y)
        return tJ2, cotJ
    end

    L2 = @view L[1:end-1, 1:end-1]
    R2 = L2'
    if Q isa StaticArray
        L2 = LowerTriangular(SMatrix{size(L2)...}(L2))
        R2 = UpperTriangular(SMatrix{size(R2)...}(R2))
    end
    y, cotJ = _find_more_nullspaces(L2, R2, y, rtol, atol, max_steps)
    tJ2 = (@view Q[:, 1:end-1]) * y

    if y isa SVector  # TODO: don't
        return SVector(tJ2...), cotJ
    end

    return tJ2, cotJ
end

function _find_more_nullspaces(L2, R2, y, rtol, atol, max_steps)
    x0 = _similar(y, length(y))
    x1 = _similar(y, length(y))

    x1 = _A_ldiv_B!(x1, R2, y)
    x1 /= norm(x1)
    y = _A_ldiv_B!(y, L2, x1)
    y /= norm(y)
    for _ in 2:max_steps
        (x0, x1) = (x1, x0)
        x1 = _A_ldiv_B!(x1, R2, y)
        x1 /= norm(x1)
        y = _A_ldiv_B!(y, L2, x1)
        y /= norm(y)
        if isapprox(x0, x1; rtol=rtol, atol=atol)
            return y, x1
        end
    end
    error("Failed to find the cokernel.")
end


function solve_simple_bifurcation!(cache, opts,
                                   u0::TV, tJ::TV,
                                   L, Q,
                                   ) where {TV <: AbstractVector}
    rtol = opts.rtol
    atol = opts.atol

    tJ1 = tJ
    tJ2, cotJ = find_more_nullspaces(Q, L, rtol, atol, opts.max_misc_steps)
    tJ2 :: TV

    function g(xi)
        prob_cache = cache.prob_cache
        v = u0 .+ xi[1] * tJ1 .+ xi[2] * tJ2
        return cotJ ⋅ residual(v, prob_cache)
    end
    hess = ForwardDiff.hessian(g, [0.0, 0.0])

    if abs(hess[1, 1]) < atol
        if abs(hess[2, 2]) < atol
            tv1 = tJ1
            tv2 = tJ2
        else
            tv1 = tJ1
            tv2 = tJ1 .* (- hess[2, 2] / 2 / hess[1, 2]) .+ tJ2
            tv2 /= norm(tv2)
        end
    else
        a = hess[1, 1]
        b = hess[1, 2]
        c = hess[2, 2]
        y1 = (- b + sqrt(b^2 - a * c)) / a
        y2 = (- b - sqrt(b^2 - a * c)) / a
        tv1 = tJ1 * y1 + tJ2
        tv2 = tJ1 * y2 + tJ2
        tv1 /= norm(tv1)
        tv2 /= norm(tv2)
    end
    @assert norm(tv1) ≈ 1
    @assert norm(tv2) ≈ 1

    if abs(tv1 ⋅ tJ1) > abs(tv2 ⋅ tJ1)
        return tv1, tv2
    else
        return tv2, tv1
    end
end


function new_branches!(cache, opts, sbint::SimpleBifurcationInterval)
    u0, tJ, L, Q = find_simple_bifurcation!(cache, opts, sbint)
    # TODO: handle not-found case

    tv1, tv2 = solve_simple_bifurcation!(cache, opts, u0, tJ, L, Q)
    tv0 = sbint.u1 .- sbint.u0
    tv0 /= norm(tv0)
    @assert abs(tv1 ⋅ tv0) > abs(tv2 ⋅ tv0)  # TODO: just warn
    tv = tv2

    args = (sbint.h, sbint.direction)
    u1, h1 = predictor_corrector_step!(cache, opts, u0, tv, args...)
    u2, h2 = predictor_corrector_step!(cache, opts, u0, -tv, args...)

    return [
        (u0, u1, sbint.direction, h1),
        (u0, u2, sbint.direction, h2),
    ]
end
