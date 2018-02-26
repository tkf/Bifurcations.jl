calc_direction(_u, J, _L, Q) = det(vcat(J, (@view Q[:, end])'))
# Q[:, end] --- tangent w/o det-fixing

find_simple_bifurcation!(cache, opts, sbint) =
    find_simple_bifurcation!(cache, opts, sbint.u0, sbint.u1,
                             sbint.direction)

find_simple_bifurcation!(cache, opts, args...) =
    find_zero!(cache, opts, calc_direction, args...)

function find_zero!(cache, opts, f, u0, u1, direction)
    prob_cache = cache.prob_cache
    H = cache.H
    J = cache.J
    Q = cache.Q
    rtol = opts.rtol
    atol = opts.atol

    H, J = residual_jacobian!(H, J, u1, prob_cache)
    A = vcat(J, _zeros(J, 1, size(J, 2)))  # TODO: improve
    L, Q = lq!(Q, A)
    f1 = f(u1, J, L, Q)

    H, J = residual_jacobian!(H, J, u0, prob_cache)
    A = vcat(J, _zeros(J, 1, size(J, 2)))  # TODO: improve
    L, Q = lq!(Q, A)
    tJ = tangent(L, Q)
    f0 = f(u0, J, L, Q)

    @assert f0 * f1 < 0

    fu = f0
    u = u0
    h = norm(u0 .- u1) / 2
    for _ in 1:opts.max_adaptations
        # predictor
        v = u .+ direction * h .* tJ

        # corrector
        for _ in 1:opts.max_corrector_steps
            w, _, H, L, Q, J = corrector_step!(H, J, Q, v, prob_cache)
            if isalmostzero(H, rtol, atol)
                cache.corrector_success = true
                @goto corrector_success
            end
            v = w
        end
        error("corrector failed")

        @label corrector_success
        tJ = tangent(L, Q)
        fv = f(v, J, L, Q)
        @assert all(isfinite.(v))

        # secant
        h = - fv / (fv - fu) * h
        u = v
        fu = fv

        if abs(h) < opts.h_zero
            return v, tJ, L, Q
        end
    end
    error("zero not found")
end
