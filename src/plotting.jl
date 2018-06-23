using RecipesBase

using .Continuations: as, ContinuationSweep,
    ContinuationSolution, sweeps_as_vectors
using .BifurcationsBase: BifurcationSweep, BifurcationSolution, special_points
using .Codim1: Codim1Sweep, Codim1Solution, Codim1Solver,
    stabilities, curves_by_stability,
    SpecialPoint, SpecialPointInterval, resolved_points

const AbstractSweep = Union{ContinuationSweep, BifurcationSweep}

domain_prototype(sweep::AbstractSweep) = as(sweep, ContinuationSweep).u[1]
domain_prototype(point::SpecialPoint) = point.u
domain_prototype(point::SpecialPointInterval) = point.u0

function var_as_index(thing, i)
    if i == 0
        i = length(domain_prototype(thing))
    end
    return i
end

const STYLE = Dict(
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
    ),
)

@recipe function plot(
        wrapper::AbstractSweep;
        vars = (0, 1),
        # To be compatible with Codim1Sweep plotter. They are ignored:
        resolve_points = nothing,
        include_points = nothing,
        )
    sweep = as(wrapper, ContinuationSweep)
    ix, iy = (var_as_index(sweep, v) for v in vars)
    xs = [u[ix] for u in sweep.u]
    ys = [u[iy] for u in sweep.u]
    (xs, ys)
end

@recipe function f(sol::ContinuationSolution; vars = (0, 1))
    ix, iy = vars
    xs = sweeps_as_vectors(sol, ix)
    ys = sweeps_as_vectors(sol, iy)

    delete!(plotattributes, :vars)
    (xs, ys)
end

@recipe function f(point::Union{SpecialPoint,
                                SpecialPointInterval};
                   vars = (0, 1),
                   style = STYLE)
    ix, iy = (var_as_index(point, v) for v in vars)

    delete!(plotattributes, :vars)
    delete!(plotattributes, :style)

    merge!(plotattributes,
           style[:point][:default],
           style[:point][point.point_type])

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

@recipe function f(sweep::Codim1Sweep;
                   vars = (0, 1),
                   resolve_points = false,
                   include_points = false)
    warn_include_points(include_points)
    ix, iy = (var_as_index(sweep, v) for v in vars)
    (info, (xs, ys)) = curves_by_stability(sweep, (ix, iy))

    # stability => style
    style = Dict(
        true => Dict(:linestyle => :solid, :linecolor => 1),
        false => Dict(:linestyle => :dot, :linecolor => 2),
    )
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

@recipe function f(sol::BifurcationSolution)
    for sweep in sol.sweeps
        @series begin
            sweep
        end
    end
end

@recipe function f(solver::BifurcationSolver;
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

function plot(plottable::Union{Codim1Sweep,
                               Codim1Solution,
                               Codim1Solver};
              resolve_points = plottable isa Codim1Solver,
              include_points = true,
              kwargs...)
    plt = Main.Plots.plot()
    for point in maybe_get_points(plottable, include_points, resolve_points)
        Main.Plots.plot!(plt, point)
    end
    Main.Plots.plot!(plt, plottable;
                     include_points = false,
                     kwargs...)
    return plt
end

plot(args...; kwargs...) = Main.Plots.plot(args...; kwargs...)
