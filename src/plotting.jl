using RecipesBase

using .Continuations: as, ContinuationSweep, ContinuationSolution,
    ContinuationSolver, sweeps_as_vectors
using .BifurcationsBase: BifurcationSweep, BifurcationSolution,
    BifurcationSolver, special_points
using .Codim1: Codim1Sweep, Codim1Solution, Codim1Solver,
    stabilities, curves_by_stability,
    SpecialPoint, SpecialPointInterval, resolved_points
using .Codim2: Codim2Sweep, Codim2Solution, Codim2Solver

const AbstractSolver = Union{ContinuationSolver, BifurcationSolver}
const AbstractSolution = Union{ContinuationSolution, BifurcationSolution}
const AbstractSweep = Union{ContinuationSweep, BifurcationSweep}
const Codim1Ctx = Union{Codim1Sweep, Codim1Solution, Codim1Solver}
const Codim2Ctx = Union{Codim2Sweep, Codim2Solution, Codim2Solver}

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

@recipe function plot(sol::BifurcationSolution)
    for sweep in sol.sweeps
        @series begin
            sweep
        end
    end
end

@recipe function plot(
        solver::BifurcationSolver;
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
            sweep
        end
    end
end

"""
    plot(...)

A thin wrapper of `Plots.plot` with some workarounds for
`Bifurcations`-related objects.
"""
function plot(plottable::Union{BifurcationSweep,
                               BifurcationSolution,
                               BifurcationSolver};
              resolve_points = plottable isa Codim1Solver,
              include_points = true,
              vars = default_vars(plottable),
              bif_style = STYLE,
              kwargs...)
    plt = Main.Plots.plot()
    for point in maybe_get_points(plottable, include_points, resolve_points)
        Main.Plots.plot!(plt, point;
                         vars = vars,
                         bif_style = bif_style)
    end
    Main.Plots.plot!(plt, plottable;
                     include_points = false,
                     vars = vars,
                     bif_style = bif_style,
                     kwargs...)
    return plt
end

plot(args...; kwargs...) = Main.Plots.plot(args...; kwargs...)
