struct NotAlmostZeroException <: Exception
    msg::AbstractString
    value::Union{AbstractArray, Number}
    rtol::Real
    atol::Real
end

function Base.showerror(io::IO, e::NotAlmostZeroException)
    n = @sprintf "%.4g" norm(e.value)
    print(e.msg, " :",
          " norm = ", n,
          " rtol = ", e.rtol,
          " atol = ", e.atol)
    if ! get(io, :compact, false)
        println(io)
        print(io, "value = ", e.value)
    end
end

function find_errors(sweep::ContinuationSweep,
                     prob::AbstractContinuationProblem,
                     opts::ContinuationOptions;
                     sweep_index = '?')
    if isempty(sweep.u)
        return []
    end
    cache = get_prob_cache(prob)
    H = _similar(sweep.u[1], length(sweep.u[1]) - 1)
    errors = []
    for (i, u) in enumerate(sweep.u)
        H = residual!(H, u, cache)
        if ! isalmostzero(H, opts.atol)
            push!(errors, NotAlmostZeroException(
                "H(sweeps[$sweep_index].u[$i]) not almost zero",
                H,
                opts.rtol,
                opts.atol,
            ))
        end
    end
    return errors
end

function find_errors(sol::ContinuationSolution, args...)
    errors = []
    for (k, sweep) in enumerate(sol.sweeps)
        append!(errors, find_errors(sweep, args...; sweep_index=k))
    end
    return errors
end

function find_errors(solver::ContinuationSolver)
    errors = find_errors(solver.sol, solver.prob, solver.opts)

    # Create a problem cache since solver.cache.prob_cache must not be
    # modified.  `find_errors` has to be a pure function!
    prob_cache = get_prob_cache(solver.prob)

    for (k, sweep) in enumerate(solver.sol.sweeps)
        if all(isindomain(u, prob_cache) for u in sweep.u)
            push!(errors, ErrorException(
                "$k-th sweep terminates before touching the boundary.",
            ))
        end
        # TODO: allow unbounded domain
    end
    return errors
end

print_errors(errors) = print_errors(stdout, errors)

function print_errors(io, errors)
    for e in errors
        showerror(io, e)
        println(io)
    end
end

# TODO: define `diagnosis` function for user-friendly API
