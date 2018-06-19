using StaticArrays: SMatrix

using ..Continuations: AbstractContinuationProblem, AbstractContinuationSolver,
    as, SweepSetup, ContinuationSweep, ContinuationSolution,
    ContinuationCache, ContinuationOptions, ContinuationSolver
import ..Continuations: step!, new_sweep!

module PointTypes
@enum(
    PointType,
    none,
    # initial,
    # simple_bifurcation,
    saddle_node,
    hopf,
    period_doubling,
)
end  # module
using .PointTypes: PointType

abstract type AbstractSpecialPoint{tkind} end

#=
struct SpecialPoint{tkind, uType, JType} <: AbstractSpecialPoint{tkind}
    point_type::PointType
    point_index::Int
    u::uType
    J::JType
    sweeps::Vector{WeakRef}
end
=#

struct SpecialPointInterval{tkind, uType, JType} <: AbstractSpecialPoint{tkind}
    timekind::tkind
    point_type::PointType
    point_index::Int
    u0::uType
    u1::uType
    # J0::JType
    J1::JType
    sweep::WeakRef
end

struct Codim1Sweep{tkind <: TimeKind,
                   S <: ContinuationSweep,
                   JType <: AbstractArray,
                   eType <: AbstractArray,
                   pType <: SpecialPointInterval{tkind}}
    super::S
    jacobians::Vector{JType}
    eigvals::Vector{eType}
    special_points::Vector{pType}
end
# See: [[../continuations/solution.jl::ContinuationSweep]]

Codim1Sweep{tkind, S, JType, eType, pType}(super::S) where{
    tkind <: TimeKind,
    S <: ContinuationSweep,
    JType <: AbstractArray,
    eType <: AbstractArray,
    pType <: SpecialPointInterval{tkind},
} = Codim1Sweep{tkind, S, JType, eType, pType}(
    super,
    JType[],
    eType[],
    pType[],
)

timekind(::Codim1Sweep{tkind}) where tkind = tkind()

eigvals_prototpye(cache::ContinuationCache) = cache.u[1:end - 1]
# TODO: improve it for SVector

function sweeptype(solver::ContinuationSolver,
                   JType::Type = typeof(solver.cache.J),
                   eType::Type = typeof(eigvals_prototpye(solver.cache)),
                   )
    tkind = typeof(timekind(solver.cache))
    S = eltype(solver.sol.sweeps)
    pType = SpecialPointInterval{tkind, eltype(S), JType}
    return Codim1Sweep{tkind, S, JType, eType, pType}
end

Base.length(sweep::Codim1Sweep) = length(as(sweep, ContinuationSweep))

struct Codim1Solution{S <: ContinuationSolution,
                      W <: Codim1Sweep}
    super::S
    sweeps::Vector{W}
end

Codim1Solution(super::ContinuationSolution,
               SweepType::Type{<: Codim1Sweep}) =
    Codim1Solution(super, SweepType[])


mutable struct Codim1Cache{P, C <: ContinuationCache{P},
                           JType, eType}
    # TODO: declare types
    super::C
    J::JType
    eigvals::eType
    point_type::PointType
end

function Codim1Cache(super::C,
                     J::JType,
                     eigvals::eType,
                     point_type = PointTypes.none,
                     ) where {P, C <: ContinuationCache{P},
                              JType, eType}
    return Codim1Cache{P, C, JType, eType}(super, J, eigvals, point_type)
end

Codim1Cache(super::ContinuationCache) =
    Codim1Cache(
        super,
        ds_jacobian(super),
        copy(eigvals_prototpye(super)),
    )

timekind(cache::Codim1Cache) = timekind(as(cache, ContinuationCache))
timekind(cache::ContinuationCache) = timekind(cache.prob_cache)

struct Codim1Solver{R <: ContinuationSolver,
                    P <: AbstractContinuationProblem,
                    C <: Codim1Cache,
                    S <: Codim1Solution,
                    } <: AbstractContinuationSolver
    super::R
    prob::P
    opts::ContinuationOptions
    cache::C
    sol::S
end

function Codim1Solver(prob::AbstractContinuationProblem,
                      opts::ContinuationOptions)
    super = ContinuationSolver(prob, opts)
    cache = Codim1Cache(super.cache)
    sol = Codim1Solution(super.sol, sweeptype(super))
    return Codim1Solver(super, prob, opts, cache, sol)
end

function push_special_point!(sweep::Codim1Sweep,
                             cache::Codim1Cache)
    push_special_point!(sweep,
                        cache.point_type,
                        as(cache, ContinuationCache).J)
end

function push_special_point!(sweep::Codim1Sweep, point_type::PointType,
                             J1)
    super = as(sweep, ContinuationSweep)
    point = SpecialPointInterval(
        timekind(sweep),
        point_type,
        length(sweep),
        super.u[end - 1],
        super.u[end],
        J1,
        WeakRef(sweep),
    )
    push!(sweep.special_points, point)
end

function new_sweep!(solver::Codim1Solver, setup::SweepSetup)
    new_sweep!(solver.super, setup)
    # calling [[./continuations/interface.jl::new_sweep!]]

    _new_sweep!(solver.sol, as(solver, ContinuationSolver))

    #=
    re_analyze!(solver, setup.u0)
    for u in setup.past_points
        re_analyze!(solver, u)
    end
    =#
end

function _new_sweep!(sol::Codim1Solution, solver::ContinuationSolver)
    super = as(sol, ContinuationSolution)
    sweep = Codim1Sweep(super.sweeps[end], solver)
    push!(sol.sweeps, sweep)
end

function Codim1Sweep(super::ContinuationSweep, solver::ContinuationSolver)
    return sweeptype(solver)(super)
end

function step!(solver::Codim1Solver)
    step!(solver.super)
    # calling [[./continuations/interface.jl::step!]]

    analyze!(solver.cache, solver.opts)
    record!(solver.sol, solver.cache)
end

function analyze!(cache::Codim1Cache, opts)
    cache.J = J = ds_jacobian(cache)
    eigvals = ds_eigvals(timekind(cache), J)
    cache.point_type = guess_point_type(timekind(cache), cache, eigvals, opts)
    cache.eigvals = eigvals
end

function record!(sol::Codim1Solution, cache::Codim1Cache)
    super = as(cache, ContinuationCache)
    sweep = sol.sweeps[end]
    push!(sweep.jacobians, copy(super.J))
    push!(sweep.eigvals, copy(cache.eigvals))
    # @assert length(sweep) == length(sweep.jacobians) == length(sweep.eigvals)

    if cache.point_type != PointTypes.none
        push_special_point!(sweep, cache)
    end
end

ds_jacobian(solver) = ds_jacobian(as(solver, ContinuationSolver).cache)
ds_jacobian(cache::Codim1Cache) = ds_jacobian(as(cache, ContinuationCache))
ds_jacobian(cache::ContinuationCache) = ds_jacobian(cache.J)
ds_jacobian(HJ::AbstractArray) = @view HJ[:, 1:end-1]
ds_jacobian(HJ::SMatrix) = HJ[:, 1:end-1]
# TOOD: optimize it for StaticArrays using generated functions
