using ..Continuations: AbstractContinuationProblem

abstract type BifurcationProblem{
        skind <: StateKind,
        tkind <: TimeKind,
    } <: AbstractContinuationProblem
end

StateKind(::Type{<: BifurcationProblem{skind}}) where skind = skind()
TimeKind(::Type{<: BifurcationProblem{_, tkind}}) where {_, tkind} = tkind()
