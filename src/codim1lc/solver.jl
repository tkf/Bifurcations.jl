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

module PointTypes
@enum(
    PointType,
    none,
)
end  # module
using .PointTypes: PointType

const Codim1LCSweep{tkind <: TimeKind, ptType <: PointType} =
    BifurcationSweep{tkind, ptType}
const Codim1LCSolution{S <: ContinuationSolution, W <: Codim1LCSweep} =
    BifurcationSolution{S, W}
# The constraints `<: TimeKind` and `<: ContinuationSolution` are
# important.
# See: [[../codim1/solver.jl::Codim1Sweep]]

BifurcationsBase.point_type_type(::Codim1LCProblem) = PointType

BifurcationsBase.regular_point(::Type{PointType}) = PointTypes.none

BifurcationsBase.eigvals_prototype(prob::Codim1LCProblem,
                                   cache::ContinuationCache) =
    cache.u[1:end - 1] .+ 1im
# TODO: improve it for SVector


mutable struct Codim1LCCache{P, C <: ContinuationCache{P},
                             JType, eType,
                             } <: BifurcationCache{P}
    super::C
    J::JType
    eigvals::eType
    point_type::PointType
end

function Codim1LCCache(super::C,
                       J::JType,
                       eigvals::eType,
                       point_type = PointTypes.none,
                       ) where {P, C <: ContinuationCache{P},
                                JType, eType}
    return Codim1LCCache{P, C, JType, eType}(super, J, eigvals, point_type)
end
# TODO: Remove this constructor after removing the type parameter `P`.

Codim1LCCache(prob::Codim1LCProblem, super::ContinuationCache) =
    Codim1LCCache(
        super,
        let v = eigvals_prototype(prob, super)
            v * v'
        end,
        # ds_jacobian(super),
        copy(eigvals_prototype(prob, super)),
    )

BifurcationsBase.BifurcationCache(prob::Codim1LCProblem,
                                  super::ContinuationCache) =
    Codim1LCCache(prob, super)

const Codim1LCSolver{
        R <: ContinuationSolver,
        P <: Codim1LCProblem,
        C <: Codim1LCCache,
        S <: Codim1LCSolution,
        } =
    BifurcationSolver{R, P, C, S}

function re_analyze!(solver::Codim1LCSolver, u::AbstractVector)
    residual_jacobian!(as(solver.cache, ContinuationCache), u)
    analyze!(solver.cache, solver.opts)

    # Suppress special point recording:
    # It's a bit ugly hack... (communicate by sharing!)
    solver.cache.point_type = PointTypes.none  # TODO: FIX!

    record!(solver.sol, solver.cache)
end

function analyze!(cache::Codim1LCCache, opts)
    #=
    cache.J = J = ds_jacobian(cache)
    eigvals = ds_eigvals(timekind(cache), J)
    cache.point_type = guess_point_type(timekind(cache), cache, eigvals, opts)
    cache.eigvals = eigvals
    =#
    set_reference!(cache)
end

function set_reference!(wrapper::Codim1LCCache)
    cache = as(wrapper, ContinuationCache)
    set_reference!(cache.prob_cache, cache.u)
end
