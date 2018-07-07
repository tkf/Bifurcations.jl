"""
    point_type_type(::BifurcationProblem) :: Type
"""
function point_type_type end  # TODO: rename

"""
    eigvals_prototype(::BifurcationProblem, ::ContinuationCache)
"""
function eigvals_prototype end

"""
    regular_point(::Type{P <: AbstractSpecialPoint}) :: P
"""
function regular_point end

"""
    analyze!(cache::BifurcationCache, opts)
"""
function analyze! end

"""
    re_analyze!(solver::BifurcationSolver, u::AbstractVector)
"""
function re_analyze! end

#=
"""
    curves(ctx [, by]) :: Vector{<: Curve}

Return a list of curves in `ctx` (sweep, solution, or solver)
optionally grouped (chunked) by a key `by` (e.g., `:stability`).
"""
function curves end
=#

"""
    measure(curve, key)
"""
function measure end
# Better names?  I thought maybe `values` but it is used for
# associative collections.  How about `observe`?  `gauge`?
# `evaluate`?  `extract`?  Maybe just define `Base.getindex`?

#=
"""
    measure(curve, keys::Tuple{T₁, ..., Tₙ}) :: Tuple{V₁, ..., Vₙ}
"""
measure(c::Any, keys::Tuple) = measure.((c,), keys)

# This makes plotting functions easy but I'm not sure if this is fine,
# especially if I were to expose it as an API.
measure(::Any, ::Nothing) = nothing
=#
