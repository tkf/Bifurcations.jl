using StaticArrays: SMatrix, SVector

function tangent(L, Q)
    tJ = Q[:, end]
    if det(Q) * det(@view L[1:end-1, 1:end-1]) < 0
        tJ *= -1
    end
    return tJ
end

function corrector_step!(H, J, Q, v, prob_cache)
    H, J = residual_jacobian!(H, J, v, prob_cache)
    A = vcat(J, _zeros(J, 1, size(J, 2)))  # TODO: improve
    L, Q = lq!(Q, A)
    y = _A_ldiv_B!((@view L[1:end-1, 1:end-1]), H)
    dv = (@view Q[:, 1:end-1]) * y
    w = v - dv
    return (w, dv, L, Q)
end

function _step!(cache, opts)
    prob_cache = cache.prob_cache
    direction = cache.direction
    u = cache.u
    H = cache.H
    J = cache.J
    Q = cache.Q
    h = cache.h
    rtol = opts.rtol
    atol = opts.atol

    cache.corrector_success = false
    cache.adaptation_success = false
    cache.simple_bifurcation = false

    H, J = residual_jacobian!(H, J, u, prob_cache)
    A = vcat(J, _zeros(J, 1, size(J, 2)))  # TODO: improve
    L, Q = lq!(Q, A)
    tJ = tangent(L, Q)

    for _ in 1:opts.max_adaptations
        # predictor
        v = u .+ direction * h .* tJ

        # corrector
        v, dv, L, Q = corrector_step!(H, J, Q, v, prob_cache)
        n1 = norm(dv)
        tJv = tangent(L, Q)
        angle = acos(min(abs(tJ ⋅ tJv), 1))  # TODO: should I use min?

        v, dv, L, Q = corrector_step!(H, J, Q, v, prob_cache)
        n2 = norm(dv)

        # step adaptation
        f_contraction =
            sqrt(n2 / n1 / opts.nominal_contraction)   # √ κ(u,h) / κ̃
        f_distance = sqrt(n1 / opts.nominal_distance)  # √ δ(u,h) / δ̃
        f_angle = angle / opts.nominal_angle_rad       #   α(u,h) / α̃
        f0 = max(
            zero_if_nan(f_contraction),
            zero_if_nan(f_distance),
            zero_if_nan(f_angle),
        )
        f = max(min(f0, 2), 1/2)
        h = h / f
        if h < opts.h_min
            return
        # elseif isalmostzero(H, rtol, atol)
        #     # If close enough to the solution, let it pass?  Should I?
        elseif f0 > 2
            continue
        end

        # corrector (again)
        for _ in 3:opts.max_corrector_steps
            if isalmostzero(H, rtol, atol)
                cache.corrector_success = true
                break
            end
            v, _, L, Q = corrector_step!(H, J, Q, v, prob_cache)
        end
        if ! cache.corrector_success
            return
        end

        cache.u = v
        cache.h = h
        cache.adaptation_success = true

        tJv = tangent(L, Q)
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

function step!(solver::ContinuationSolver, max_steps)
    cache = solver.cache
    for _ in 1:max_steps
        step!(solver)
        if ! isindomain(cache.u, cache.prob_cache)
            return true
        end
    end
    return false
end
