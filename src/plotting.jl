using Statistics: mean

using Parameters: @with_kw, @unpack
using RecipesBase
using Setfield: @set, print_application, Lens, PropertyLens

using .ArrayUtils: nan_

using .Continuations: as, ContinuationSweep, ContinuationSolution,
    ContinuationSolver, sweeps_as_vectors
using .BifurcationsBase: BifurcationSweep, BifurcationSolution,
    BifurcationSolver, special_points, measure, problem_of
using .Codim1: Codim1Sweep, Codim1Solution, Codim1Solver,
    AbstractCodim1SpecialPoint, stabilities, curves_by_stability,
    SpecialPoint, SpecialPointInterval, resolved_points
using .Codim2: Codim2Sweep, Codim2Solution, Codim2Solver,
    AbstractCodim2SpecialPoint
using .Codim1LimitCycle: Codim1LCSweep, Codim1LCSolution, Codim1LCSolver,
    CWStateMeasurement, limitcycles
using .Codim2LimitCycle: Codim2LCSweep, Codim2LCSolution, Codim2LCSolver

const AbstractSolver = Union{ContinuationSolver, BifurcationSolver}
const AbstractSolution = Union{ContinuationSolution, BifurcationSolution}
const AbstractSweep = Union{ContinuationSweep, BifurcationSweep}
const Codim1Ctx = Union{Codim1Sweep, Codim1Solution, Codim1Solver,
                        AbstractCodim1SpecialPoint}
const Codim2Ctx = Union{Codim2Sweep, Codim2Solution, Codim2Solver,
                        AbstractCodim2SpecialPoint}
const Codim1LCCtx = Union{Codim1LCSweep, Codim1LCSolution, Codim1LCSolver}
const Codim2LCCtx = Union{Codim2LCSweep, Codim2LCSolution, Codim2LCSolver}
const LCCtx = Union{Codim1LCCtx, Codim2LCCtx}

@with_kw struct StateSpacePlotter{TC, TM}
    ctx::TC
    mapping::TM
end

StateSpacePlotter(ctx) =
    StateSpacePlotter(ctx, (:x => 1, :y => 2))
StateSpacePlotter(ctx::LCCtx) =
    StateSpacePlotter(ctx, (:x => 1, :y => 2, :color => :period))

plot_state_space(args...; kwargs...) =
    RecipesBase.plot(StateSpacePlotter(args...); kwargs...)

plot_state_space!(args...; kwargs...) =
    RecipesBase.plot!(StateSpacePlotter(args...); kwargs...)

plottables(ctx, key) = plottables(measure(ctx, key))
plottables(ctx, keys::Tuple) = plottables.((ctx,), keys)
plottables(::Any, ::Nothing) = nothing
# plottables(ctx, keys) = [plottables(c, keys) for c in curves(ctx)]

# post-processing:
plottables(x::Number) = x
plottables(xs::AbstractArray) = xs
plottables(xs::Vector{<: Tuple}) = hcat(collect.(xs)...)'

process_key(::Any, key::Integer) = key
mapping_to_tuple(::Any; x, y) = (x, y)

"""
    get_keys(ctx, mapping; kwargs...) :: Tuple

Transform `mapping` to measurement keys.
"""
get_keys(ctx, mapping; kwargs...) =
    process_key.((ctx,), mapping_to_tuple(ctx; mapping...);
                 kwargs...)
# TODO: Maybe get rid of get_keys?  There are too much indirections here.

map_to_plottables(ctx, mapping; kwargs...) =
    plottables(ctx, get_keys(ctx, mapping; kwargs...))

dim_domain(ctx) = length(domain_prototype(ctx))
domain_prototype(solver::AbstractSolver) = domain_prototype(solver.sol)
domain_prototype(sol::AbstractSolution) = domain_prototype(sol.sweeps[1])
domain_prototype(sweep::AbstractSweep) = as(sweep, ContinuationSweep).u[1]
domain_prototype(point::SpecialPoint) = point.u
domain_prototype(point::SpecialPointInterval) = point.u0

# TODO: don't define process_key(::Codim1Ctx, ...); use measure() interface.
process_key(::Codim1Ctx, key::Integer) = key
function process_key(ctx::Codim1Ctx, key::Symbol)
    if key in (:p1, :parameter)
        return length(domain_prototype(ctx))
    else
        error("Unsupported mapping key: ", key)
    end
end
# - [[./fixedpoint.jl::^function.*get_u0]]

# TODO: don't define process_key(::Codim2Ctx, ...); use measure() interface.
process_key(::Codim2Ctx, key::Integer) = key
function process_key(ctx::Codim2Ctx, key::Symbol)
    if key == :p1
        return length(domain_prototype(ctx)) - 1
    elseif key == :p2
        return length(domain_prototype(ctx))
    elseif key in (:ω, :w)
        if ! (contkind(problem_of(ctx)) isa HopfCont)
            error("Key $key is invalid for $(contkind(prob))")
        end
        return length(domain_prototype(ctx)) - 2
    else
        error("Unsupported mapping key: ", key)
    end
end
# - [[./codim2/diffeq.jl::^function.*get_u0]]
# - [[./codim2/diffeq.jl::get(param_axis]]

assert_non_empty(solver::AbstractSolver) = assert_non_empty(solver.sol)
function assert_non_empty(sol::AbstractSolution)
    if isempty(sol.sweeps)
        error("Solution does not have a sweep:\n",
              sprint(showcompact, sol))
    end
    for sweep in sol.sweeps
        assert_non_empty(sweep)
    end
end
function assert_non_empty(sweep::AbstractSweep)
    if isempty(as(sweep, ContinuationSweep).u)
        error("Sweep does not have a point:\n",
              sprint(showcompact, sweep))
    end
end

default_mapping(::Any) = (:x => 1, :y => 2)
default_mapping(::Codim1Ctx) = (:x => :p1, :y => 1)
default_mapping(::Codim1LCCtx) = (:x => :p1, :y => :p2, :color => :stability)
default_mapping(::Union{Codim2Ctx, Codim2LCCtx}) = (:x => :p1, :y => :p2)

lens_name(l::Lens) = sprint(print_application, l)
lens_name(l::PropertyLens) = lstrip(invoke(lens_name, Tuple{Lens}, l), '.')

function var_name(ctx::Codim1Ctx, i::Integer)
    prob = problem_of(ctx)
    if i == dim_domain(ctx)
        return lens_name(prob.p.param_axis)
    else
        return "u$i"
    end
end
# - [[./fixedpoint.jl::^function.*get_u0]]

function var_name(ctx::Codim2Ctx, i::Integer)
    prob = problem_of(ctx)
    if i == dim_domain(ctx)
        return lens_name(prob.param_axis2)
    elseif i == dim_domain(ctx) - 1
        return lens_name(prob.param_axis1)
    elseif contkind(prob) isa HopfCont && i == dim_domain(ctx) - 2
        return "ω"
    else
        return "u$i"
    end
end
# - [[./codim2/diffeq.jl::^function.*get_u0]]
# - [[./codim2/diffeq.jl::get(param_axis]]

const STYLE = Dict(
    :stability => Dict(
        # is_stable => plotattributes
        true => Dict(:linestyle => :solid, :linecolor => 1),
        false => Dict(:linestyle => :dot, :linecolor => 2),
    ),
    :point => Dict(
        :default => Dict(
            :markersize => 5,
        ),
        Codim1.PointTypes.saddle_node => Dict(
            :markercolor => 3,
            :markershape => :cross,
        ),
        Codim1.PointTypes.hopf => Dict(
            :markercolor => 4,
            :markershape => :circle,
        ),
        Codim1.PointTypes.period_doubling => Dict(
            :markercolor => 5,
            :markershape => :triangle,
        ),
        Codim2.PointTypes.cusp => Dict(
            :markercolor => 6,
            :markershape => :diamond,
        ),
        Codim2.PointTypes.bautin => Dict(
            :markercolor => 7,
            :markershape => :circle,
        ),
        Codim2.PointTypes.bogdanov_takens => Dict(
            :markercolor => 8,
            :markershape => :hexagon,
        ),
        Codim2.PointTypes.fold_hopf => Dict(
            :markercolor => 9,
            :markershape => :star8,
        ),
        Codim2.PointTypes.hopf_hopf => Dict(
            :markercolor => 10,
            :markershape => :star5,
        ),
    ),
)

struct PlottableData{T, K}
    data::T
    keys::K
end

@recipe function plot(
        pltbl::PlottableData,
        )
    color_key = pltbl.keys[3]
    if color_key !== nothing
        clims = nan_(extrema, (
            zs
            for data in pltbl.data
            for zs in plottables(data, color_key)
        ))
        clims --> clims
        colorbar_title --> string(color_key)
    end
    for data in pltbl.data
        xs, ys, zs = plottables(data, pltbl.keys)
        @series begin
            label --> ""
            if zs !== nothing
                line_z := zs
                marker_z := zs
            end
            (xs, ys)
        end
    end
end

@recipe function plot(
        wrapper::AbstractSweep,
        mapping = default_mapping(wrapper);
        # To be compatible with Codim1Sweep plotter. They are ignored:
        bif_style = nothing,
        resolve_points = nothing,
        include_points = nothing,
        )
    sweep = as(wrapper, ContinuationSweep)
    ix, iy = get_keys(sweep, mapping) :: Tuple{<: Integer, <: Integer}
    xs = [u[ix] for u in sweep.u]
    ys = [u[iy] for u in sweep.u]
    (xs, ys)
end

@recipe function plot(sol::ContinuationSolution,
                      mapping)
    ix, iy = get_keys(sol, mapping) :: Tuple{<: Integer, <: Integer}
    xs = sweeps_as_vectors(sol, ix)
    ys = sweeps_as_vectors(sol, iy)

    (xs, ys)
end

plottables(point::SpecialPoint, key::Integer) = [point.u[key]]
plottables(point::SpecialPointInterval, key::Integer) =
    [point.u0[key], point.u1[key]]
# TODO: plot SpecialPointInterval as an interval

@recipe function plot(
        point::Union{SpecialPoint,
                     SpecialPointInterval},
        mapping = default_mapping(point);
        bif_style = STYLE)

    if haskey(bif_style, :point)
        point_style = bif_style[:point]
        merge!(plotattributes,
               get(point_style, :default, Dict()),
               get(point_style, point.point_type, Dict()))
    end

    label --> ""
    map_to_plottables(point, mapping)
end

function maybe_get_points(sweep, include_points, resolve_points)
    if include_points
        if resolve_points
            return resolved_points(sweep)
        else
            return special_points(sweep)
        end
    end
    return []
end

@recipe function plot(
        sweep::Codim1Sweep,
        mapping = default_mapping(sweep);
        bif_style = STYLE,
        resolve_points = false,
        include_points = true)

    ix, iy = get_keys(sweep, mapping) :: Tuple{<: Integer, <: Integer}
    (info, (xs, ys)) = curves_by_stability(sweep, (ix, iy))

    # is_stable => Dict(:linestyle => ?, :linecolor => ?)
    style = deepcopy(STYLE[:stability])
    # Empty `style` won't work so copying the default...  Not ideal.
    # TODO: support arbitrary keys, not just linestyle and linecolor
    if haskey(bif_style, :stability)
        stability_style = bif_style[:stability]
        merge!(style[true], get(stability_style, true, Dict()))
        merge!(style[false], get(stability_style, false, Dict()))
    end

    linestyle = []
    linecolor = []
    for stable in info
        push!(linestyle, style[stable][:linestyle])
        push!(linecolor, style[stable][:linecolor])
    end

    for point in maybe_get_points(sweep, include_points, resolve_points)
        @series begin
            point
        end
    end

    xname = var_name(sweep, ix)
    yname = var_name(sweep, iy)
    for (x, y, lc, ls) in zip(xs, ys, linecolor, linestyle)
        @series begin
            label --> ""
            linecolor --> lc
            linestyle --> ls
            xlabel --> xname
            ylabel --> yname
            (x, y)
        end
    end
end

mapping_to_tuple(
    ::LCCtx;
    x,
    y,
    color = nothing,
    # color,  # TODO: use it
) = (x, y, color)

process_key(::LCCtx, key::Any; _...) = key
process_key(::LCCtx, key::Integer; cw_agg=nothing) =
    if cw_agg !== nothing
        CWStateMeasurement(key, cw_agg)
    else
        key
    end

@recipe function plot(
        sweep::Union{Codim1LCSweep, Codim2LCSweep},
        mapping = default_mapping(sweep);
        cw_agg = extrema,  # coordinate-wise aggregator
        resolve_points = false,
        include_points = false,
        bif_style = STYLE)

    # TODO: support special points
    #=
    for point in maybe_get_points(sweep, include_points, resolve_points)
        @series begin
            (point, mapping)
        end
    end
    =#

    keys = get_keys(sweep, mapping; cw_agg = cw_agg)
    curve_list = [sweep]  # TODO: chunk by stability
    # curve_list = curves(sweep, :stability)

    PlottableData(curve_list, keys)
end

@recipe function plot(
        ssp::StateSpacePlotter{<: LCCtx},
        bif_style = STYLE)
    @unpack ctx, mapping = ssp

    # TODO: support special "points"

    keys = get_keys(ctx, mapping)

    xlabel --> "u$(keys[1])"
    ylabel --> "u$(keys[2])"
    colorbar_title --> "$(keys[3])"
    PlottableData(limitcycles(ctx), keys)
end

@recipe function plot(
        sweep::Codim2Sweep,
        mapping = default_mapping(sweep);
        resolve_points = false,
        include_points = true)

    for point in maybe_get_points(sweep, include_points, resolve_points)
        @series begin
            (point, mapping)
        end
    end

    ix, iy = get_keys(sweep, mapping) :: Tuple{<: Integer, <: Integer}
    xname = var_name(sweep, ix)
    yname = var_name(sweep, iy)

    label --> ""
    linecolor --> 1
    xlabel --> xname
    ylabel --> yname
    as(sweep, ContinuationSweep), (x=ix, y=iy)  # plot(::AbstractSweep)
end

@recipe function plot(
        sol::BifurcationSolution,
        arguments__...)  # TODO: use `mapping`
    for sweep in sol.sweeps
        @series begin
            # `args` instead of `arguments__` does not work
            (sweep, arguments__...)
        end
    end
end

@recipe function plot(
        solver::BifurcationSolver,
        arguments__...;  # TODO: use `mapping`
        resolve_points = solver isa Union{Codim1Solver, Codim2Solver},
        include_points = true)

    for point in maybe_get_points(solver, include_points, resolve_points)
        @series begin
            point
        end
    end

    for sweep in solver.sol.sweeps
        @series begin
            include_points := false  # already included
            (sweep, arguments__...)
        end
    end
end

#=
@recipe function plot(
        ssp::StateSpacePlotter{<: BifurcationSolution},
        )
    for sweep in ssp.ctx.sweeps
        @series begin
            @set ssp.ctx = sweep
        end
    end
end

@recipe function plot(
        ssp::StateSpacePlotter{<: BifurcationSolver},
        )
    @set ssp.ctx = ssp.ctx.sol
end
=#
