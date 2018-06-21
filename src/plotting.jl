using RecipesBase

using .Continuations: as, ContinuationSweep,
    ContinuationSolution, sweeps_as_vectors
using .Codim1: Codim1Sweep, Codim1Solution, stabilities, curves_by_stability,
    SpecialPoint, SpecialPointInterval, special_points, resolved_points

const AbstractSweep = Union{ContinuationSweep, Codim1Sweep}

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

@recipe function f(sweep::Codim1Sweep;
                   vars = (0, 1),
                   resolve_points = false,
                   include_points = false)
    # TODO: turn on `include_points`.  It is turned off at the moment
    # since it disturbs line style and color.
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

    if include_points
        if resolve_points
            point_list = resolved_points(sweep)
        else
            point_list = special_points(sweep)
        end

        for point in point_list
            @series begin
                point
            end
        end
    end

    label --> ""
    linecolor --> reshape(linecolor, (1, length(linecolor)))
    linestyle --> reshape(linestyle, (1, length(linestyle)))
    (xs, ys)
end

@recipe function f(sol::Codim1Solution)
    for sweep in sol.sweeps
        @series begin
            sweep
        end
    end
end
