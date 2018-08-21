special_intervals(solver::BifurcationSolver, args...) =
    special_intervals(solver.sol, args...)

special_intervals(sol::BifurcationSolution, args...) =
    [point for sweep in sol.sweeps for point in special_intervals(sweep, args...)]

special_intervals(sweep::BifurcationSweep) = sweep.special_intervals
special_intervals(sweep::BifurcationSweep, point_type::Enum) =
    special_intervals(sweep, (point_type,))
special_intervals(sweep::BifurcationSweep,
               point_types::NTuple{N, <:Enum}) where N =
    filter(p -> p.point_type in point_types, special_intervals(sweep))


"""
    problem_of(ctx) :: BifurcationProblem

Get the problem associated with `ctx` (a solver, solution, or sweep).
"""
problem_of(solver::ContinuationSolver) = solver.prob
problem_of(sol::ContinuationSolution) = sol.prob
problem_of(sweep::ContinuationSweep) = problem_of(sweep.sol.value)

problem_of(solver::BifurcationSolver) = solver.prob
problem_of(sol::BifurcationSolution) = problem_of(as(sol, ContinuationSolution))
problem_of(sweep::BifurcationSweep) = problem_of(as(sweep, ContinuationSweep))

get_param_axis(ctx, ::Val{1}) = _get_param_axis1(contkind(ctx), ctx)
_get_param_axis1(::OneParamCont, ctx) = problem_of(ctx).param_axis
_get_param_axis1(::FixedPointCont, ctx) = problem_of(ctx).p.param_axis  # TODO: remove
_get_param_axis1(::TwoParamCont, ctx) = problem_of(ctx).param_axis1

get_param_axis(ctx, ::Val{2}) = problem_of(ctx).param_axis2
get_param_axis1(ctx) = get_param_axis(ctx, Val(1))
get_param_axis2(ctx) = get_param_axis(ctx, Val(2))

function count_special_points(sweep::BifurcationSweep)
    counter = Dict{point_type_type(sweep), Int}()
    for sp in sweep.special_intervals
        counter[sp.point_type] = get!(counter, sp.point_type, 0) + 1
    end
    return counter
end

function count_special_points(sol::BifurcationSolution)
    counter = Dict{PointTypeType(eltype(sol.sweeps)), Int}()
    for sweep in sol.sweeps
        for (k, v) in count_special_points(sweep)
            counter[k] = get!(counter, k, 0) + v
        end
    end
    return counter
end

count_special_points(solver::BifurcationSolver) =
    count_special_points(solver.sol)
