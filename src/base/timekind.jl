abstract type TimeKind end
struct Discrete <: TimeKind end
struct Continuous <: TimeKind end

timekind(::T) where T = TimeKind(T)


using ..Continuations: AbstractContinuationCache
TimeKind(::Type{<: AbstractContinuationCache{PC}}) where PC = TimeKind(PC)
