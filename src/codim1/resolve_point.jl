using Parameters: @unpack
using ..Continuations: find_zero!
using ..BifurcationsBase: AbstractSpecialPoint

resolved_points(solver::Codim1Solver) =
    [resolve_point(point, solver) for point in special_points(solver)]

function resolve_point(point::AbstractSpecialPoint, solver::Codim1Solver)
    super = as(solver, ContinuationSolver)
    return resolve_point(point, super.cache, super.opts)
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
    f = functional_for(point_type, timekind(cache))
    direction = as(point.sweep.value, ContinuationSweep).direction
    u, tJ, L, Q, J =
        find_zero!(cache, opts, f, point.u0, point.u1, direction)
    return SpecialPoint(timekind(point),
                        point_type, point_index, u, J,
                        WeakRef(point.sweep.value))
end

function functional_for(point_type::PointType, tkind::TimeKind)
    ptype = Val{point_type}()
    # For call signature of `f`, see:
    # [[../continuations/zero_point.jl::f(]]
    return (u, J, L, Q) -> functional_for(ptype, tkind, u, J, L, Q)
end
# It evaluates to a scalar so I'm calling it a functional.  Are there
# better names?  Root-findee?

const Instability = Union{
    Val{PointTypes.saddle_node},
    Val{PointTypes.period_doubling},
    Val{PointTypes.hopf},
}

functional_for(::Instability, tkind::TimeKind, u, J, L, Q) =
    stability_index(tkind, ds_jacobian(J))

stability_index(tkind::Discrete, J) = ds_eigvals(tkind, J)[1] - 1
stability_index(tkind::Continuous, J) = ds_eigvals(tkind, J)[1] # TODO: improve!
# I don't need to calculate all eigenvalues here.
