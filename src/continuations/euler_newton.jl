using StaticArrays: SMatrix, SVector

# TODO: improve tangent
tangent(A::Matrix) = nullspace(A)[:, 1]
tangent(A::SMatrix) = SVector(tangent(Matrix(A))...)

# TODO: improve _pinv
_pinv(A::Matrix) = pinv(A)
_pinv(A::SMatrix{S1, S2}) where {S1, S2} = SMatrix{S2, S1}(pinv(Matrix(A)))

function corrector_step!(H, J, v, prob_cache)
    H, J = residual_jacobian!(H, J, v, prob_cache)
    dv = _pinv(J) * H
    w = v - dv
    return (H, J, w, dv)
end

function _step!(cache, opts)
    prob_cache = cache.prob_cache
    direction = cache.direction
    u = cache.u
    H = cache.H
    J = cache.J
    h = cache.h
    rtol = opts.rtol
    atol = opts.atol

    cache.corrector_success = false
    cache.adaptation_success = false
    cache.simple_bifurcation = false

    H, J = residual_jacobian!(H, J, u, prob_cache)
    tJ = tangent(J)

    for _ in 1:opts.max_adaptations
        # predictor
        v = u .+ direction * h .* tJ

        # corrector
        H, J, v, dv = corrector_step!(H, J, v, prob_cache)
        n1 = norm(dv)
        tJv = tangent(J)
        angle = acos(min(abs(tJ ⋅ tJv), 1))  # TODO: should I use min?

        H, J, v, dv = corrector_step!(H, J, v, prob_cache)
        n2 = norm(dv)

        # step adaptation
        f0 = max(
            sqrt(n2 / n1 / opts.nominal_contraction),  # √ κ(u,h) / κ̃
            sqrt(n1 / opts.nominal_distance),          # √ δ(u,h) / δ̃
            angle / opts.nominal_angle_rad,            #   α(u,h) / α̃
        )
        f = max(min(f0, 2), 1/2)
        h = h / f
        if h < opts.h_min
            return
        elseif f0 > 2
            continue
        end

        # corrector (again)
        for _ in 3:opts.max_corrector_steps
            if isalmostzero(H, rtol, atol)
                cache.corrector_success = true
                break
            end
            H, J, v, _ = corrector_step!(H, J, v, prob_cache)
        end
        if ! cache.corrector_success
            return
        end

        cache.u = v
        cache.h = h
        cache.adaptation_success = true

        tJv = tangent(J)
        if tJ ⋅ tJv < 0
            cache.direction *= -1
            cache.simple_bifurcation = true
        end
        return
    end
end

function record!(sol, cache)
    push_point!(sol, cache.u, cache.simple_bifurcation)
end

function step!(solver::ContinuationSolver)
    _step!(solver.cache, solver.opts)
    if ! solver.cache.adaptation_success
        error("Failed to adapt steplength h.")
    end
    if ! solver.cache.corrector_success
        error("Failed in corrector loop.")
    end
    record!(solver.sol, solver.cache)
    solver.i += 1
end
