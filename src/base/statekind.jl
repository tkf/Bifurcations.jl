abstract type StateKind end
struct MutableState <: StateKind end
struct ImmutableState <: StateKind end

statekind(::T) where T = StateKind(T)
# TODO: Use StateKind everywhere instead of boolean iip type parameter.


StateKind(::Type{<: AbstractContinuationCache{PC}}) where PC = StateKind(PC)
StateKind(::Type{<: AbstractProblemCache{P}}) where P = StateKind(P)
