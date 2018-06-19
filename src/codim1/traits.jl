abstract type TimeKind end
struct Discrete <: TimeKind end
struct Continuous <: TimeKind end

timekind(::T) where T = TimeKind(T)
