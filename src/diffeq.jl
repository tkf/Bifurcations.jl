using DiffEqBase: AbstractODEProblem, DiscreteProblem
using Setfield: Lens, set, get

const DEP{iip} = AbstractODEProblem{uType, tType, iip} where {uType, tType}


abstract type StateKind end
struct MutableState <: StateKind end
struct ImmutableState <: StateKind end

statekind(::T) where T = StateKind(T)
# TODO: Move StateKind to continuations/base.jl and use it everywhere.
# It's only used in here at the moment.


TimeKind(::Type{<: DiscreteProblem}) = Discrete()
TimeKind(::Type{<: AbstractODEProblem}) = Continuous()
StateKind(::Type{<: DEP{true}}) = MutableState()
StateKind(::Type{<: DEP{false}}) = ImmutableState()

struct DiffEqWrapper{P, L}
    de_prob::P
    param_axis::L
end

function diffeq_homotopy(H, x, p::DiffEqWrapper{<:DEP{true}}, t)
    q = set(p.param_axis, p.de_prob.p, t)
    p.de_prob.f(H, x, q, 0)
    maybe_subtract!(H, x, statekind(p.de_prob), timekind(p.de_prob))
end

function diffeq_homotopy(x, p::DiffEqWrapper{<:DEP{false}}, t)
    q = set(p.param_axis, p.de_prob.p, t)
    H = p.de_prob.f(x, q, 0)
    return maybe_subtract!(H, x, statekind(p.de_prob), timekind(p.de_prob))
end

maybe_subtract!(H, ::Any, ::StateKind, ::Continuous) = H
maybe_subtract!(H, x, ::MutableState, ::Discrete) = H .-= x
maybe_subtract!(H, x, ::ImmutableState, ::Discrete) = H .- x


"""
    FixedPointBifurcationProblem(ode_or_map::AbstractODEProblem,
                                 param_axis::Lens,
                                 t_domain::Tuple;
                                 <keyword arguments>)

# Arguments
- `ode_or_map`: An `ODEProblem` or `DiscreteProblem`.
- `param_axis :: Lens`: The lens to set/get a parameter of `ode_or_map`.
- `t_domain :: Tuple`: A pair of numbers specifying the lower and
  upper bound for `param_axis`.
"""
function FixedPointBifurcationProblem(
        de_prob::DEP{iip}, param_axis::Lens, t_domain::Tuple;
        kwargs...) where iip
    u0 = de_prob.u0
    t0 = get(param_axis, de_prob.p)
    p = DiffEqWrapper(deepcopy(de_prob), param_axis)
    return FixedPointBifurcationProblem{iip, typeof(timekind(de_prob))}(
        diffeq_homotopy, u0, t0,
        t_domain, p; kwargs...)
end
