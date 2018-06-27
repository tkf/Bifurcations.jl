using Parameters: @unpack
using ..Continuations: find_zero!
using ..BifurcationsBase: AbstractSpecialPoint, special_points,
    contkind, FixedPointCont, ContinuationKind

resolved_points(solver::BifurcationSolver) =
    [resolve_point(point, solver) for point in special_points(solver)]

resolved_points(solver::BifurcationSolver,
                point_types::NTuple{N, <:Enum},
                ) where {N} =
    [resolve_point(point, solver) for point in special_points(solver)
     if point.point_type in point_types]

resolved_points(solver::BifurcationSolver, point_type::Enum) =
    resolved_points(solver, (point_type,))

function resolve_point(point::AbstractSpecialPoint, solver::BifurcationSolver)
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
    f = testfn_for(point_type, timekind(cache), contkind(cache))
    direction = as(point.sweep.value, ContinuationSweep).direction
    u, tJ, L, Q, J =
        find_zero!(cache, opts, f, point.u0, point.u1, direction)
    return SpecialPoint(timekind(point),
                        point_type, point_index, u, J,
                        WeakRef(point.sweep.value))
end

function testfn_for(point_type::Enum, tkind::TimeKind, ckind::ContinuationKind)
    ptype = Val{point_type}()
    # For call signature of `f`, see:
    # [[../continuations/zero_point.jl::f(]]
    return (u, J, L, Q) -> testfn(ptype, tkind, ckind, u, J, L, Q)
end

const Instability = Union{
    Val{PointTypes.saddle_node},
    Val{PointTypes.period_doubling},
    Val{PointTypes.hopf},
}

testfn(::Instability, tkind::TimeKind, ::FixedPointCont,
       u, J, L, Q) =
    stability_index(tkind, ds_jacobian(J))

stability_index(tkind::Discrete, J) =
    abs(ds_eigvals(tkind, J)[1]) - 1
stability_index(tkind::Continuous, J) =
    real(ds_eigvals(tkind, J)[1]) # TODO: improve!
# I don't need to calculate all eigenvalues here.
