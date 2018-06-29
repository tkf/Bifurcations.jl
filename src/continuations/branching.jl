using ForwardDiff
using StaticArrays: SMatrix, SArray
using Setfield: @set


"""
Find vectors tJ2 and cotJ such that ker(J) = span{tJ1, tJ2} and
ker(J') = span{cotJ} (left nullspace or cokernel):
"""
function find_more_nullspaces(Q, L, rtol, atol, max_steps)
    y = _zeros(Q, size(Q, 1) - 1)
    if y isa SArray
        y = @set y[length(y)] = 1   # y[end] doesn't work
    else
        y[end] = 1
    end

    #=
    if abs(L[end-1, end-1]) < atol
        tJ2 = Q[:, end-1]
        cotJ = (@view Q[1:end-1, 1:end-1]) \ y
        return tJ2, cotJ
    end
    =#

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

function _find_more_nullspaces(L2, R2,
                               y :: T,
                               rtol, atol, max_steps,
                               ) where {T <: AbstractVector}
    if abs(det(L2)) < atol
        # TODO: optimize
        ker_L2 = nullspace(Array(L2))
        ker_R2 = nullspace(Array(R2))
        if size(ker_L2, 2) > 0 && size(ker_R2, 2) > 0
            return (T(@view ker_L2[:, 1]), T(@view ker_R2[:, 1]))
        end
        # Otherwise, let's fallback to the manual method.
    end
    x0 = _similar(y, length(y))
    x1 = _similar(y, length(y))

    x1 = _normalize!(R2 \ y)
    y = _normalize!(L2 \ x1)
    for _ in 2:max_steps
        (x0, x1) = (x1, x0)
        x1 = _normalize!(R2 \ y)
        y = _normalize!(L2 \ x1)
        if isapprox(x0, x1; rtol=rtol, atol=atol)
            return y::T, x1::T
        end
    end
    error("Failed to find the cokernel.")
end


@with_kw struct SimpleBifurcationSolution{RV, LV, M}
    found::Bool
    tv1::RV
    tv2::RV
    tJ1::RV
    tJ2::RV
    cotJ::LV
    hess::M
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
    H11, H12, _, H22 = hess

    # Solve (ξ₁, ξ₂)ᵀ H (ξ₁, ξ₂) = 0 where ξ₁=x and ξ₂=1-x
   discriminant = - H11 * H22 + H12^2  # = det(hess)
    if discriminant <= 0
        tv1 = tJ1 .* NaN
        tv2 = tJ2 .* NaN
        return SimpleBifurcationSolution(false, tv1, tv2, tJ1, tJ2, cotJ, hess)
    end

    x1 = (- H12 + H22 + sqrt(discriminant)) / (H11 + H22 - 2 * H12)
    x2 = (- H12 + H22 - sqrt(discriminant)) / (H11 + H22 - 2 * H12)
    if abs(x1) < abs(x2)
        x1, x2 = x2, x1
    end
    tv1 = @. x1 * tJ1 + (1 - x1) * tJ2
    tv2 = @. x2 * tJ1 + (1 - x2) * tJ2
    tv1 /= norm(tv1)
    tv2 /= norm(tv2)

    return SimpleBifurcationSolution(true, tv1, tv2, tJ1, tJ2, cotJ, hess)
end


function new_branches!(cache, opts, sbint::SimpleBifurcationInterval)
    u0, tJ, L, Q = find_simple_bifurcation!(cache, opts, sbint)
    # TODO: handle not-found case

    sbsol = solve_simple_bifurcation!(cache, opts, u0, tJ, L, Q)
    @unpack tv1, tv2 = sbsol
    @assert sbsol.found  # TODO: nicer error handling

    # Choose the direction `tv` of the new branch.  Use the one least
    # parallel to the direction `tv0` along the previous curve.
    tv0 = sbint.u1 .- sbint.u0
    tv0 /= norm(tv0)
    if abs(tv1 ⋅ tv0) > abs(tv2 ⋅ tv0)
        tv = tv2
    else
        tv = tv1
    end

    args = (sbint.h, sbint.direction)
    u1, h1 = predictor_corrector_step!(cache, opts, u0, tv, args...)
    u2, h2 = predictor_corrector_step!(cache, opts, u0, -tv, args...)

    # TODO: Make error message creation lazy; i.e., store the points
    # and vectors in a custom exception type and construct the message
    # in showerror.
    parallel_to_old = du -> abs(du ⋅ tv0) > abs(du ⋅ tv)
    failures = []
    if parallel_to_old(u1 .- u0) || parallel_to_old(u2 .- u0)
        push!(failures, string(
            "New branch candidates are more parallel to the old curve",
            " than the computed new direction."))
    end
    if (isapprox(u1, u0; atol=opts.atol, rtol=opts.rtol) ||
        isapprox(u2, u0; atol=opts.atol, rtol=opts.rtol))
        push!(failures,
              "New points are approximately equal to the branch point.")
    end
    if ! isempty(failures)
        failure_msg = join(failures, "\n    * ")
        warn("""
            Failed to branch off from old curve.
            Reason(s):
                * $(failure_msg)
            """)
        # Not sure about this anymore:
        #=
            Possible fix(es):
                * Decrease h0 (= $(opts.h0))
        =#
    end

    return [
        (u0, u1, sbint.direction, h1),
        (u0, u2, sbint.direction, h2),
    ]
end
