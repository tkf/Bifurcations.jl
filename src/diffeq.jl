using DiffEqBase: AbstractODEProblem, DiscreteProblem
using Setfield: Lens, set, get, @set

const DEP{iip} = AbstractODEProblem{uType, tType, iip} where {uType, tType}


TimeKind(::Type{<: DiscreteProblem}) = Discrete()
TimeKind(::Type{<: AbstractODEProblem}) = Continuous()
StateKind(::Type{<: DEP{true}}) = MutableState()
StateKind(::Type{<: DEP{false}}) = ImmutableState()

struct DiffEqWrapper{P, L, D}
    de_prob::P
    param_axis::L
    param_diff::D
end

DiffEqWrapper(de_prob, param_axis) = DiffEqWrapper(de_prob, param_axis, nothing)

function diffeq_homotopy(H, x, p::DiffEqWrapper{<:DEP{true}}, t)
    q = set(p.de_prob.p, p.param_axis, t)
    p.de_prob.f(H, x, q, 0)
    maybe_subtract!(H, x, statekind(p.de_prob), timekind(p.de_prob))
end

function diffeq_homotopy(x, p::DiffEqWrapper{<:DEP{false}}, t)
    q = set(p.de_prob.p, p.param_axis, t)
    H = p.de_prob.f(x, q, 0)
    return maybe_subtract!(H, x, statekind(p.de_prob), timekind(p.de_prob))
end

maybe_subtract!(H, x, p) = maybe_subtract!(H, x, statekind(p), timekind(p))
maybe_subtract!(H, ::Any, ::StateKind, ::Continuous) = H
maybe_subtract!(H, x, ::MutableState, ::Discrete) = H .-= x
maybe_subtract!(H, x, ::ImmutableState, ::Discrete) = H .- x

struct MutableParamDiff{F, X, C}
    param_to_state::F
    x::X
    cfg::C
end

function MutableParamDiff(p::DiffEqWrapper)
    x = deepcopy(p.de_prob.u0)
    param_to_state = let x = x
        (H, t) -> diffeq_homotopy(H, x, p, t)
    end
    cfg = ForwardDiff.DerivativeConfig(param_to_state, x, 1.0)
    return MutableParamDiff(param_to_state, x, cfg)
end

function diffeq_homotopy_jacobian(H, J, x, p::DiffEqWrapper{<:DEP{true}}, t)
    p.param_diff :: MutableParamDiff
    p.param_diff.x .= x
    ForwardDiff.derivative!(
        (@view J[:, end]),
        p.param_diff.param_to_state,
        H,
        t,
        p.param_diff.cfg,
    )
    q = set(p.de_prob.p, p.param_axis, t)
    p.de_prob.f.jac((@view J[:, 1:end-1]), x, q, t)
    return (H, maybe_subtract!(J, Eye(size(J)...), p.de_prob))
end

"""
    BifurcationProblem(ode_or_map::AbstractODEProblem,
                       param_axis::Lens,
                       t_domain::Tuple;
                       <keyword arguments>)

# Arguments
- `ode_or_map`: An `ODEProblem` or `DiscreteProblem`.
- `param_axis :: Lens`: The lens to set/get a parameter of `ode_or_map`.
- `t_domain :: Tuple`: A pair of numbers specifying the lower and
  upper bound for `param_axis`.
"""
function BifurcationProblem(
        de_prob::DEP, param_axis::Lens, t_domain::Tuple;
        kwargs0...)
    de_prob = deepcopy(de_prob)
    u0 = de_prob.u0
    t0 = get(de_prob.p, param_axis)
    p0 = DiffEqWrapper(de_prob, param_axis)
    if de_prob.f.jac === nothing
        p = p0
        kwargs = kwargs0
    elseif de_prob isa DEP{true}
        p = @set p0.param_diff = MutableParamDiff(p0)
        kwargs = (; homotopy_jacobian=diffeq_homotopy_jacobian, kwargs0...)
    else
        # TODO: implement
        p = p0
        kwargs = kwargs0
    end
    return FixedPointBifurcationProblem(
        statekind(de_prob),
        timekind(de_prob),
        diffeq_homotopy, u0, t0,
        t_domain, p; kwargs...)
end
