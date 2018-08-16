abstract type ContinuationKind end
abstract type OneParamCont <: ContinuationKind end
abstract type TwoParamCont <: ContinuationKind end

# Codim1
struct FixedPointCont <: OneParamCont end

# Codim2
struct SaddleNodeCont <: TwoParamCont end
struct HopfCont <: TwoParamCont end

# Codim1LimitCycle
struct LimitCycleCont <: OneParamCont end

# Codim2LimitCycle
struct FoldLimitCycleCont <: TwoParamCont end


#=
contkind(::T) where T = ContinuationKind(T)
=#

function contkind(cache::AbstractContinuationCache)  # TODO: trait
    prob = as(cache, ContinuationCache).prob_cache.prob  # TODO: interface
    return contkind(prob)
end

contkind(solver::AbstractContinuationSolver) = contkind(solver.prob)
