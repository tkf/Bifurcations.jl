using RecipesBase

using .Continuations: as, ContinuationSweep,
    ContinuationSolution, sweeps_as_vectors

function var_as_index(wrapper, i)
    sweep = as(wrapper, ContinuationSweep)
    if i == 0
        i = length(sweep.u[1])
    end
    return i
end

@recipe function f(sol::ContinuationSolution; vars = (0, 1))
    ix, iy = vars
    xs = sweeps_as_vectors(sol, ix)
    ys = sweeps_as_vectors(sol, iy)

    delete!(plotattributes, :vars)
    (xs, ys)
end

using .Codim1: Codim1Sweep, Codim1Solution, stabilities, curves_by_stability

@recipe function f(sweep::Codim1Sweep; vars = (0, 1))
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
