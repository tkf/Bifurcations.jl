using Parameters: @with_kw, @unpack
using RecipesBase
using Setfield: @set

using .CompatUtils: @required
using .ArrayUtils: nan_

using .Continuations: as, ContinuationSweep, ContinuationSolution,
    ContinuationSolver, sweeps_as_vectors
using .BifurcationsBase: BifurcationSweep, BifurcationSolution,
    BifurcationSolver, special_points, measure
using .Codim1: Codim1Sweep, Codim1Solution, Codim1Solver,
    stabilities, curves_by_stability,
    SpecialPoint, SpecialPointInterval, resolved_points
using .Codim2: Codim2Sweep, Codim2Solution, Codim2Solver
using .Codim1LimitCycle: Codim1LCSweep, Codim1LCSolution, Codim1LCSolver,
    CWStateMeasurement, limitcycles

const AbstractSolver = Union{ContinuationSolver, BifurcationSolver}
const AbstractSolution = Union{ContinuationSolution, BifurcationSolution}
const AbstractSweep = Union{ContinuationSweep, BifurcationSweep}
const Codim1Ctx = Union{Codim1Sweep, Codim1Solution, Codim1Solver}
const Codim2Ctx = Union{Codim2Sweep, Codim2Solution, Codim2Solver}
const Codim1LCCtx = Union{Codim1LCSweep, Codim1LCSolution, Codim1LCSolver}

@with_kw struct StateSpacePlotter{TC, TM}
    ctx::TC
    mapping::TM
end

StateSpacePlotter(ctx) =
    StateSpacePlotter(ctx, (:x => 1, :y => 2))
StateSpacePlotter(ctx::Codim1LCCtx) =
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

get_keys(ctx, mapping) = _get_keys(ctx; mapping...)

dim_domain(ctx) = length(domain_prototype(ctx))
domain_prototype(solver::AbstractSolver) = domain_prototype(solver.sol)
domain_prototype(sol::AbstractSolution) = domain_prototype(sol.sweeps[1])
domain_prototype(sweep::AbstractSweep) = as(sweep, ContinuationSweep).u[1]
domain_prototype(point::SpecialPoint) = point.u
domain_prototype(point::SpecialPointInterval) = point.u0

function var_as_index(thing, i)
    if i == 0
        i = length(domain_prototype(thing))
    end
    return i
end

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

default_vars(::Any) = (0, 1)
function default_vars(ctx::Codim1Ctx)
    assert_non_empty(ctx)
    D = dim_domain(ctx)
    return (D, 1)
end
function default_vars(ctx::Codim2Ctx)
    assert_non_empty(ctx)
    D = dim_domain(ctx)
    return (D - 1, D)
end
default_vars(ctx::Codim1LCCtx) = (:x => :p1, :y => 1)

function _merge_to_mapping(ctx, mapping, vars)   # TODO: remove
    # Get rid of `vars` in all plotting functions.
    if vars !== nothing
        @assert mapping === nothing
        return vars
    elseif mapping === nothing
        return default_vars(ctx)
    end
    return mapping
end

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
        wrapper::AbstractSweep;
        vars = default_vars(wrapper),
        # To be compatible with Codim1Sweep plotter. They are ignored:
        bif_style = nothing,
        resolve_points = nothing,
        include_points = nothing,
        )
    sweep = as(wrapper, ContinuationSweep)
    ix, iy = (var_as_index(sweep, v) for v in vars)
    xs = [u[ix] for u in sweep.u]
    ys = [u[iy] for u in sweep.u]
    (xs, ys)
end

@recipe function plot(sol::ContinuationSolution;
                      vars = default_vars(sol))
    ix, iy = vars
    xs = sweeps_as_vectors(sol, ix)
    ys = sweeps_as_vectors(sol, iy)

    delete!(plotattributes, :vars)
    (xs, ys)
end

@recipe function plot(
        point::Union{SpecialPoint,
                     SpecialPointInterval};
        vars = default_vars(point),
        bif_style = STYLE)
    ix, iy = (var_as_index(point, v) for v in vars)

    if haskey(bif_style, :point)
        point_style = bif_style[:point]
        merge!(plotattributes,
               get(point_style, :default, Dict()),
               get(point_style, point.point_type, Dict()))
    end

    label --> ""
    if point isa SpecialPointInterval
        xs = [point.u0[ix], point.u1[ix]]
        ys = [point.u0[iy], point.u1[iy]]
        # TODO: plot interval
        ([mean(xs)], [mean(ys)])
    else
        ([point.u[ix]], [point.u[iy]])
    end
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

function warn_include_points(include_points)
    if include_points
        warn("""include_points = true is set.
             Note that it is known to disturb line styles and colors
             such that stability information is wrongly plotted.
             As a workaround, use:
                 Bifurcations.plot(sweep)
                 Bifurcations.plot(sol)
                 Bifurcations.plot(solver)
             """)
    end
end
# TODO: Make `include_points` work.

@recipe function plot(
        sweep::Codim1Sweep;
        vars = default_vars(sweep),
        bif_style = STYLE,
        resolve_points = false,
        include_points = false)
    warn_include_points(include_points)
    ix, iy = (var_as_index(sweep, v) for v in vars)
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

    label --> ""
    linecolor --> reshape(linecolor, (1, length(linecolor)))
    linestyle --> reshape(linestyle, (1, length(linestyle)))
    (xs, ys)
end

_get_keys(
    ::Codim1LCCtx;
    (@required x),
    (@required y),
    color = nothing,
    # (@required color),  # TODO: use it
) = (x, y, color)

process_key(::Codim1LCCtx, key::Any; _...) = key
process_key(::Codim1LCCtx, key::Integer; cw_agg=nothing) =
    if cw_agg !== nothing
        CWStateMeasurement(key, cw_agg)
    else
        key
    end

@recipe function plot(
        sweep::Codim1LCSweep,
        mapping = nothing;
        # mapping = (:x => :p1, :y => 1, :color => :stability);  # TODO: use it
        vars = nothing,
        cw_agg = extrema,  # coordinate-wise aggregator
        resolve_points = false,
        include_points = false,
        bif_style = STYLE)

    mapping = _merge_to_mapping(sweep, mapping, vars)

    warn_include_points(include_points)
    # TODO: support special points; note that I have to get rid of
    # `vars` first.
    #=
    for point in maybe_get_points(sweep, include_points, resolve_points)
        @series begin
            (point, mapping)
        end
    end
    =#

    keys = process_key.((sweep,), get_keys(sweep, mapping);
                        cw_agg = cw_agg)

    curve_list = [sweep]  # TODO: chunk by stability
    # curve_list = curves(sweep, :stability)

    PlottableData(curve_list, keys)
end

@recipe function plot(
        ssp::StateSpacePlotter{<: Codim1LCCtx},
        vars = nothing,
        bif_style = STYLE)
    @unpack ctx, mapping = ssp
    mapping = _merge_to_mapping(ctx, mapping, vars)

    # TODO: support special "points"

    keys = get_keys(ctx, mapping)

    PlottableData(limitcycles(ctx), keys)
end

@recipe function plot(
        sweep::Codim2Sweep;
        resolve_points = false,
        include_points = false)

    for point in maybe_get_points(sweep, include_points, resolve_points)
        @series begin
            point
        end
    end

    label --> ""
    linecolor --> 1
    as(sweep, ContinuationSweep)  # plot(::AbstractSweep)
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
        resolve_points = solver isa Codim1Solver,
        include_points = true)
    warn_include_points(include_points)

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

"""
    plot!(...)

A thin wrapper of `Plots.plot!` with some workarounds for
`Bifurcations`-related objects.
"""
function plot!(plt,
               plottable::Union{BifurcationSweep,
                                BifurcationSolution,
                                BifurcationSolver},
               args...;
               resolve_points = plottable isa Codim1Solver,
               include_points = true,
               vars = length(args) > 0 ? nothing : default_vars(plottable),
               bif_style = STYLE,
               kwargs...)
    for point in maybe_get_points(plottable, include_points, resolve_points)
        Main.Plots.plot!(plt, point, args...;
                         vars = vars,
                         bif_style = bif_style)
    end
    Main.Plots.plot!(plt, plottable, args...;
                     include_points = false,
                     vars = vars,
                     bif_style = bif_style,
                     kwargs...)
    return plt
end

plot!(args...; kwargs...) = Main.Plots.plot!(args...; kwargs...)


"""
    plot(...)

A thin wrapper of `Plots.plot` with some workarounds for
`Bifurcations`-related objects.
"""
function plot(plottable::Union{BifurcationSweep,
                               BifurcationSolution,
                               BifurcationSolver},
              args...;
              kwargs...)
    plt = Main.Plots.plot()
    plot!(plt, plottable, args...; kwargs...)
    return plt
end

plot(args...; kwargs...) = Main.Plots.plot(args...; kwargs...)
