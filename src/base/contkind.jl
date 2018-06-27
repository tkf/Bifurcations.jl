abstract type ContinuationKind end

# Codim1
struct FixedPointCont <: ContinuationKind end

# Codim2
struct SaddleNodeCont <: ContinuationKind end
struct HopfCont <: ContinuationKind end


#=
contkind(::T) where T = ContinuationKind(T)
=#

function contkind(cache::AbstractContinuationCache)  # TODO: trait
    prob = as(cache, ContinuationCache).prob_cache.prob  # TODO: interface
    return contkind(prob)
end
