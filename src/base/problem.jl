using ..Continuations: AbstractContinuationProblem

abstract type BifurcationProblem{
        skind <: StateKind,
        tkind <: TimeKind,
    } <: AbstractContinuationProblem
end
