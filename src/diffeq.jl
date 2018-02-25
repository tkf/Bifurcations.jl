using DiffEqBase: AbstractODEProblem
using Setfield: Lens, set, get

const DEP{iip} = AbstractODEProblem{uType, tType, iip} where {uType, tType}

struct DiffEqWrapper{P, L}
    de_prob::P
    param_axis::L
end

function diffeq_homotopy(H, x, p::DiffEqWrapper{<:DEP{true}}, t)
    q = set(p.param_axis, p.de_prob.p, t)
    p.de_prob.f(H, x, q, 0)
end

function diffeq_homotopy(x, p::DiffEqWrapper{<:DEP{false}}, t)
    q = set(p.param_axis, p.de_prob.p, t)
    return p.de_prob.f(x, q, 0)
end

function FixedPointBifurcationProblem(
        de_prob::DEP{iip}, param_axis::Lens, t_domain::Tuple;
        kwargs...) where iip
    u0 = de_prob.u0
    t0 = get(param_axis, de_prob.p)
    p = DiffEqWrapper(deepcopy(de_prob), param_axis)
    return FixedPointBifurcationProblem{iip}(diffeq_homotopy, u0, t0,
                                             t_domain, p; kwargs...)
end
