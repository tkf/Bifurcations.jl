using RecipesBase

using .Continuations: as, ContinuationSolution, sweeps_as_vectors

@recipe function f(sol::ContinuationSolution; vars = (0, 1))
    ix, iy = vars
    xs = sweeps_as_vectors(sol, ix)
    ys = sweeps_as_vectors(sol, iy)

    delete!(plotattributes, :vars)
    (xs, ys)
end

using .Codim1: Codim1Solution

@recipe function f(sol::Codim1Solution; vars = (0, 1))
    super = as(sol, ContinuationSolution)
    super
end
