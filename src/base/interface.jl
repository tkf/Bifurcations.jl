"""
    point_type_type(::BifurcationProblem) :: Type
"""
function point_type_type end  # TODO: rename

"""
    eigvals_prototpye(::BifurcationProblem, ::ContinuationCache)
"""
function eigvals_prototpye end

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
