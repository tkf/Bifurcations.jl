using Parameters: @unpack
using ..Continuations: find_zero!
using ..BifurcationsBase: AbstractSpecialPoint, special_intervals,
    special_points, push_special_point!,
    contkind, FixedPointCont, ContinuationKind

function default_resolve_exception_handler(err, interval, warn_exceptions)
    if any(isa.(Ref(err), warn_exceptions))
        @warn(
            """
            Failed to find bifurcation point within:
            $(interval)
            """,
            exception = err)
        return true
    end
    return false
end

default_resolve_exception_handler(warn_exceptions::Tuple) =
    (args...) -> default_resolve_exception_handler(args..., warn_exceptions)

function resolving_points(
        solver::BifurcationSolver, args...;
        warn_exceptions = (SingularException,),
        exception_handler = default_resolve_exception_handler(warn_exceptions),
        )
    Iterators.filter(
        point -> point !== nothing,
        try
            resolve_point(interval, solver)
        catch err
            # [[../continuations/branching.jl::SingularException]]
            if ! exception_handler(err, interval)
                rethrow()
            end
            nothing
        end for interval in special_intervals(solver, args...))
end

resolved_points(solver::BifurcationSolver, args...; kwargs...) =
    collect(resolving_points(solver, args...; kwargs...))

function resolve_point(point::AbstractSpecialPoint, solver::BifurcationSolver)
    super = as(solver, ContinuationSolver)
    return resolve_point(point, super.cache, super.opts)
end

"""
    resolve_points!(solver; kwargs...)

Resolve special points if not.  It is a no-op if all points are
already resolved.
"""
function resolve_points!(
        solver::BifurcationSolver;
        warn_exceptions = (SingularException,),
        exception_handler = default_resolve_exception_handler(warn_exceptions),
        )
    for sweep in solver.sol.sweeps
        ofs = length(sweep.special_points) + 1
        for interval in @view sweep.special_intervals[ofs:end]
            point = try
                resolve_point!(interval,
                               as(solver.cache, ContinuationCache),
                               solver.opts)
            catch err
                if !exception_handler(err, interval)
                    rethrow()
                end
                continue
            end
            push_special_point!(sweep, point)
        end
    end
end

function special_points!(solver, args...; kwargs...)
    resolve_points!(solver; kwargs...)
    return special_points(solver, args...)
end

resolve_point(point::AbstractSpecialPoint,
              cache::ContinuationCache,
              args...) =
    resolve_point!(point, deepcopy(cache), args...)

resolve_point!(point::SpecialPoint, _...) = point  # already resolved

function resolve_point!(point::SpecialPointInterval,
                        cache::ContinuationCache,
                        opts::ContinuationOptions,
                        ) :: SpecialPoint
    @unpack point_type, point_index = point
    f = testfn_for(point_type, timekind(cache), contkind(cache),
                   cache.prob_cache)
    direction = as(point.sweep.value, ContinuationSweep).direction
    u, tJ, L, Q, J =
        find_zero!(cache, opts, f, point.u0, point.u1, direction)
    return SpecialPoint(timekind(point),
                        point_type, point_index, u, J,
                        WeakRef(point.sweep.value))
end

function testfn_for(point_type::Enum, tkind::TimeKind, ckind::ContinuationKind,
                    prob_cache)
    ptype = Val{point_type}()
    # For call signature of `f`, see:
    # [[../continuations/zero_point.jl::f(]]
    return (u, J, L, Q) -> testfn(ptype, tkind, ckind, prob_cache, u, J, L, Q)
end

testfn(::Val{PointTypes.saddle_node}, ::Continuous, ::FixedPointCont,
       prob_cache, u, J, L, Q) =
    det(ds_jacobian(J))

function testfn(::Val{PointTypes.hopf}, ::Continuous, ::FixedPointCont,
                prob_cache, u, J, L, Q)
    vals = _eigvals(ds_jacobian(J))
    _, i = findmin(abs.(real.(vals)))
    return real(vals[i])
end

testfn(::Val{PointTypes.saddle_node}, ::Discrete, ::FixedPointCont,
       prob_cache, u, J, L, Q) =
    det(ds_jacobian(J) - I)

testfn(::Val{PointTypes.period_doubling}, ::Discrete, ::FixedPointCont,
       prob_cache, u, J, L, Q) =
    det(ds_jacobian(J) + I)
