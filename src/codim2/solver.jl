using StaticArrays: SMatrix

using ..Continuations: AbstractContinuationProblem, AbstractContinuationSolver,
    AbstractContinuationCache,
    as, SweepSetup, ContinuationSweep, ContinuationSolution,
    ContinuationCache, ContinuationOptions, ContinuationSolver,
    residual_jacobian!
import ..Continuations: step!, new_sweep!

using ..BifurcationsBase: timekind, Continuous, Discrete
import ..BifurcationsBase: TimeKind

using ..BifurcationsBase: contkind, SaddleNodeCont, HopfCont

using ..BifurcationsBase
using ..BifurcationsBase: SpecialPoint, SpecialPointInterval,
    BifurcationSweep, BifurcationSolution, BifurcationSolver,
    BifurcationProblem, BifurcationCache,
    eigvals_prototype, allocate_sweep!, check_sweep_length, record!
import ..BifurcationsBase: analyze!, re_analyze!


"""
Codimension 2 bifurcation points and special points.

These points are supported:

- [Cusp bifurcation]
- [Bogdanov-Takens bifurcation]
- ([Bautin Bifurcation] not yet)
- [Fold-Hopf Bifurcation]
- [Hopf-Hopf Bifurcation]

[Cusp bifurcation]:
http://www.scholarpedia.org/article/Cusp_bifurcation
[Bogdanov-Takens bifurcation]:
http://www.scholarpedia.org/article/Bogdanov-Takens_bifurcation
[Bautin Bifurcation]:
http://www.scholarpedia.org/article/Bautin_bifurcation
[Fold-Hopf Bifurcation]:
http://www.scholarpedia.org/article/Fold-Hopf_bifurcation
[Hopf-Hopf Bifurcation]:
http://www.scholarpedia.org/article/Hopf-Hopf_bifurcation

* http://www.scholarpedia.org/article/Bifurcation
* https://www.encyclopediaofmath.org/index.php/Codimension-two_bifurcations
"""
module PointTypes
@enum(
    PointType,
    none,
    # initial,
    cusp,
    bautin,
    bogdanov_takens,
    fold_hopf,
    hopf_hopf,
)
end  # module
using .PointTypes: PointType

const Codim2Sweep{tkind <: TimeKind, ptType <: PointType} =
    BifurcationSweep{tkind, ptType}
const Codim2Solution{S <: ContinuationSolution, W <: Codim2Sweep} =
    BifurcationSolution{S, W}
# The constraints `<: TimeKind` and `<: ContinuationSolution` are
# important.
# See: [[../codim1/solver.jl::Codim1Sweep]]

BifurcationsBase.point_type_type(::Codim2Problem) = PointType

BifurcationsBase.regular_point(::Type{PointType}) = PointTypes.none

function BifurcationsBase.eigvals_prototype(prob::Codim2Problem,
                                            cache::ContinuationCache)
    N = length(cache.u) ÷ 2 - 1
    return cache.u[1:N] .+ 1im
end
# TODO: improve it for SVector


mutable struct Codim2Cache{P, C <: ContinuationCache{P},
                           JType, eType,
                           } <: BifurcationCache{P}
    super::C
    J::JType
    eigvals::eType
    prev_eigvals::eType
    point_type::PointType

    """
    The quadratic coefficient ``a(0)``.
    See: http://www.scholarpedia.org/article/Saddle-node_bifurcation
    """
    quadratic_coefficient::Float64
    prev_quadratic_coefficient::Float64
end
# TODO: Maybe rename Codim2Cache to SNContCache for something.
# Or put bifurcation type-specific cache in a sub-cache?

function Codim2Cache(super::C,
                     J::JType,
                     eigvals::eType,
                     point_type::PointType = PointTypes.none,
                     ) where {P, C <: ContinuationCache{P},
                              JType, eType}
    return Codim2Cache{P, C, JType, eType}(super, J,
                                           eigvals, copy(eigvals),
                                           point_type,
                                           NaN, NaN)
end
# TODO: Remove this constructor after removing the type parameter `P`.

Codim2Cache(prob::Codim2Problem, super::ContinuationCache) =
    Codim2Cache(
        super,
        ds_jacobian(super),
        copy(eigvals_prototype(prob, super)),
    )

BifurcationsBase.BifurcationCache(prob::Codim2Problem,
                                  super::ContinuationCache) =
    Codim2Cache(prob, super)

const Codim2Solver{
        R <: ContinuationSolver,
        P <: Codim2Problem,
        C <: Codim2Cache,
        S <: Codim2Solution,
        } =
    BifurcationSolver{R, P, C, S}

# TODO: trait
function BifurcationsBase.contkind(prob::DiffEqCodim2Problem)
    if eltype(prob.v0) <: Complex
        return HopfCont()
    else
        return SaddleNodeCont()
    end
end

function re_analyze!(solver::Codim2Solver, u::AbstractVector)
    residual_jacobian!(as(solver.cache, ContinuationCache), u)
    analyze!(solver.cache, solver.opts)

    # Suppress special point recording:
    # It's a bit ugly hack... (communicate by sharing!)
    solver.cache.point_type = PointTypes.none  # TODO: FIX!

    record!(solver.sol, solver.cache)
end

function analyze!(cache::Codim2Cache, opts)
    cnt = contkind(cache)
    cache.J = J = ds_jacobian(cache)
    cache.prev_eigvals = cache.eigvals
    cache.eigvals = ds_eigvals(timekind(cache), J)
    set_quadratic_coefficient!(cnt, cache)
    cache.point_type = guess_point_type(cnt, timekind(cache), cache, opts)
    set_augsys_cache!(cache)
end

ds_jacobian(solver) = ds_jacobian(as(solver, ContinuationSolver).cache)
ds_jacobian(cache::Codim2Cache) = ds_jacobian(as(cache, ContinuationCache))

function ds_jacobian(cache::ContinuationCache)
    prob = cache.prob_cache.prob  # TODO: interface
    return ds_jacobian(prob, cache.J)
end

ds_jacobian(prob::DiffEqCodim2Problem, J) =
    Array(_ds_jacobian(eltype(prob.v0), J))
# Workaround: Only hermitian matrices are diagonalizable by *StaticArrays*.

function _ds_jacobian(E::Type, J::AbstractMatrix)
    d = dims_from_augsys(size(J, 2), E)
    return @view J[1:d.ds_dim, 1:d.ds_dim]
end

@generated function _ds_jacobian(::Type{E},
                                 J::SMatrix{S1, S2, T}
                                 ) where {E, S1, S2, T}
    d = dims_from_augsys(S2, E)
    N = d.ds_dim
    values = [:(J[$i, $j]) for j in 1:N for i in 1:N]
    quote
        return SMatrix{$N, $N, $T}($(values...))
    end
end

function ds_state(cache::ContinuationCache)
    prob = cache.prob_cache.prob  # TODO: interface
    return ds_state(prob, cache.u)
end

# TODO: define ds_f!
ds_f(x, prob, cache::Codim2Cache) = ds_f(x, as(cache, ContinuationCache))

function ds_f(x, cache::ContinuationCache)
    prob = cache.prob_cache.prob  # TODO: interface
    f = prob.de_prob.f  # TODO: interface
    p = modified_param!(prob, cache.u)
    return f(x, p, 0)
end

set_quadratic_coefficient!(::HopfCont, ::Codim2Cache) = nothing

function set_quadratic_coefficient!(::SaddleNodeCont, cache::Codim2Cache)
    cache.prev_quadratic_coefficient = cache.quadratic_coefficient
    cache.quadratic_coefficient =
        sn_quadratic_coefficient(timekind(cache), statekind(cache), cache)
end
