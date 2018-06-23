using StaticArrays: SMatrix

using ..Continuations: AbstractContinuationProblem, AbstractContinuationSolver,
    AbstractContinuationCache,
    as, SweepSetup, ContinuationSweep, ContinuationSolution,
    ContinuationCache, ContinuationOptions, ContinuationSolver,
    residual_jacobian!
import ..Continuations: step!, new_sweep!

using ..BifurcationsBase: timekind, Continuous, Discrete
import ..BifurcationsBase: TimeKind

using ..BifurcationsBase
using ..BifurcationsBase: SpecialPoint, SpecialPointInterval,
    BifurcationSweep, BifurcationSolution, BifurcationSolver,
    BifurcationProblem, BifurcationCache,
    eigvals_prototype, allocate_sweep!, check_sweep_length, record!
import ..BifurcationsBase: analyze!, re_analyze!


"""
Codimension 2 bifurcation points and special points.

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
    return cache.u[1:N]
end
# TODO: improve it for SVector


mutable struct Codim2Cache{P, C <: ContinuationCache{P},
                           JType, eType,
                           } <: BifurcationCache{P}
    super::C
    J::JType
    eigvals::eType
    point_type::PointType
end

function Codim2Cache(super::C,
                     J::JType,
                     eigvals::eType,
                     point_type = PointTypes.none,
                     ) where {P, C <: ContinuationCache{P},
                              JType, eType}
    return Codim2Cache{P, C, JType, eType}(super, J, eigvals, point_type)
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

function re_analyze!(solver::Codim2Solver, u::AbstractVector)
    residual_jacobian!(as(solver.cache, ContinuationCache), u)
    analyze!(solver.cache, solver.opts)

    # Suppress special point recording:
    # It's a bit ugly hack... (communicate by sharing!)
    solver.cache.point_type = PointTypes.none  # TODO: FIX!

    record!(solver.sol, solver.cache)
end

function analyze!(cache::Codim2Cache, opts)
    cache.J = J = ds_jacobian(cache)
    cache.point_type = guess_point_type(timekind(cache), cache, opts)
    cache.eigvals = ds_eigvals(timekind(cache), J)
end

ds_jacobian(solver) = ds_jacobian(as(solver, ContinuationSolver).cache)
ds_jacobian(cache::Codim2Cache) = ds_jacobian(as(cache, ContinuationCache))
ds_jacobian(cache::ContinuationCache) = ds_jacobian(cache.J)

function ds_jacobian(HJ::AbstractArray)
    N = size(HJ, 2) ÷ 2 - 1
    return @view HJ[1:N, 1:N]
end

function ds_jacobian(HJ::SMatrix)
    N = size(HJ, 2) ÷ 2 - 1
    return HJ[1:N, 1:N]
end
# TOOD: optimize it for StaticArrays using generated functions
