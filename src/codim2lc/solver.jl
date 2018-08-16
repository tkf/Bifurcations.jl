using ..Continuations: AbstractContinuationProblem, AbstractContinuationSolver,
    AbstractContinuationCache,
    as, SweepSetup, ContinuationSweep, ContinuationSolution,
    ContinuationCache, ContinuationOptions, ContinuationSolver,
    residual_jacobian!

using ..BifurcationsBase: timekind, Continuous, Discrete
import ..BifurcationsBase: TimeKind

using ..BifurcationsBase
using ..BifurcationsBase: SpecialPoint, SpecialPointInterval,
    BifurcationSweep, BifurcationSolution, BifurcationSolver,
    BifurcationProblem, BifurcationCache,
    eigvals_prototype, allocate_sweep!, check_sweep_length, record!
import ..BifurcationsBase: analyze!, re_analyze!

using ..Codim1LimitCycle: set_reference!

module PointTypes
@enum(
    PointType,
    none,
)
end  # module
using .PointTypes: PointType

const Codim2LCSweep{tkind <: TimeKind, ptType <: PointType} =
    BifurcationSweep{tkind, ptType}
const Codim2LCSolution{S <: ContinuationSolution, W <: Codim2LCSweep} =
    BifurcationSolution{S, W}
# The constraints `<: TimeKind` and `<: ContinuationSolution` are
# important.
# See: [[../codim1/solver.jl::Codim1Sweep]]

BifurcationsBase.point_type_type(::Codim2LCProblem) = PointType

BifurcationsBase.regular_point(::Type{PointType}) = PointTypes.none

function BifurcationsBase.eigvals_prototype(prob::Codim2LCProblem,
                                            cache::ContinuationCache)
    n = length(prob.super.de_prob.u0)
    return similar(cache.u[1:1] .+ 1im, n)
end


mutable struct Codim2LCCache{P, C <: ContinuationCache{P},
                             JType, eType,
                             } <: BifurcationCache{P}
    super::C
    J::JType
    eigvals::eType
    point_type::PointType
end

function Codim2LCCache(super::C,
                       J::JType,
                       eigvals::eType,
                       point_type = PointTypes.none,
                       ) where {P, C <: ContinuationCache{P},
                                JType, eType}
    return Codim2LCCache{P, C, JType, eType}(super, J, eigvals, point_type)
end
# TODO: Remove this constructor after removing the type parameter `P`.

Codim2LCCache(prob::Codim2LCProblem, super::ContinuationCache) =
    Codim2LCCache(
        super,
        let v = eigvals_prototype(prob, super)
            v * v'
        end,
        # ds_jacobian(super),
        copy(eigvals_prototype(prob, super)),
    )

BifurcationsBase.BifurcationCache(prob::Codim2LCProblem,
                                  super::ContinuationCache) =
    Codim2LCCache(prob, super)

const Codim2LCSolver{
        R <: ContinuationSolver,
        P <: Codim2LCProblem,
        C <: Codim2LCCache,
        S <: Codim2LCSolution,
        } =
    BifurcationSolver{R, P, C, S}

# TODO: trait
BifurcationsBase.contkind(::Codim2LCProblem) = FoldLimitCycleCont()

function re_analyze!(solver::Codim2LCSolver, u::AbstractVector)
    residual_jacobian!(as(solver.cache, ContinuationCache), u)
    analyze!(solver.cache, solver.opts)

    # Suppress special point recording:
    # It's a bit ugly hack... (communicate by sharing!)
    solver.cache.point_type = PointTypes.none  # TODO: FIX!

    record!(solver.sol, solver.cache)
end

function analyze!(cache::Codim2LCCache, opts)
    #=
    cache.J = J = ds_jacobian(cache)
    eigvals = ds_eigvals(timekind(cache), J)
    cache.point_type = guess_point_type(timekind(cache), cache, eigvals, opts)
    cache.eigvals = eigvals
    =#

    # This works since the "layout" of the head of cache.super.u is
    # identical to `LimitCycleCache` case.  Maybe better to not rely
    # on the layout implicitly:
    set_reference!(cache)
end
