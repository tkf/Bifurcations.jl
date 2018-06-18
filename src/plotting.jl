using RecipesBase

@recipe function f(sol::ContinuationSolution; vars = (0, 1))
    ix, iy = vars
    xs = sweeps_as_vectors(sol, ix)
    ys = sweeps_as_vectors(sol, iy)

    delete!(plotattributes, :vars)
    (xs, ys)
end
