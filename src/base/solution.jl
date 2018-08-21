using ..Continuations: as,
    ContinuationSweep, ContinuationSolution,
    ContinuationCache, ContinuationSolver

abstract type AbstractSpecialPoint{tkind <: TimeKind, ptType} end
TimeKind(::Type{<: AbstractSpecialPoint{tkind}}) where tkind = tkind()

struct SpecialPoint{tkind, ptType, uType, JType,
                    } <: AbstractSpecialPoint{tkind, ptType}
    timekind::tkind
    point_type::ptType
    point_index::Int
    u::uType
    J::JType
    sweep::WeakRef
end
# TODO: put BifurcationProblem in SpecialPoint/Interval struct so that
# init(point) can initiate codim+1 bifurcation solver.

struct SpecialPointInterval{tkind, ptType, uType, JType,
                            } <: AbstractSpecialPoint{tkind, ptType}
    timekind::tkind
    point_type::ptType
    point_index::Int
    u0::uType
    u1::uType
    # J0::JType
    J1::JType
    sweep::WeakRef
end

struct BifurcationSweep{
        tkind <: TimeKind,
        ptType,
        S <: ContinuationSweep,
        JType <: AbstractArray,
        eType <: AbstractArray,
        siType <: SpecialPointInterval{tkind, ptType},
        spType <: SpecialPoint{tkind, ptType},
        }
    timekind::tkind
    super::S
    jacobians::Vector{JType}
    eigvals::Vector{eType}
    special_intervals::Vector{siType}
    special_points::Vector{spType}
end
# See: [[../continuations/solution.jl::ContinuationSweep]]

BifurcationSweep{tkind, ptType, S, JType, eType, siType, spType}(super::S) where{
    tkind <: TimeKind,
    ptType,
    S <: ContinuationSweep,
    JType <: AbstractArray,
    eType <: AbstractArray,
    siType <: SpecialPointInterval{tkind, ptType},
    spType <: SpecialPoint{tkind, ptType},
} = BifurcationSweep{tkind, ptType, S, JType, eType, siType, spType}(
    tkind(),
    super,
    JType[],
    eType[],
    siType[],
    spType[],
)

TimeKind(::Type{<: BifurcationSweep{tkind}}) where tkind = tkind()

Base.length(sweep::BifurcationSweep) = length(as(sweep, ContinuationSweep))

PointTypeType(::Type{<: BifurcationSweep{_tk, ptType}}) where {_tk, ptType} =
    ptType

point_type_type(::T) where {T <: BifurcationSweep} = PointTypeType(T)

contkind(sweep::BifurcationSweep) = contkind(problem_of(sweep)) # TODO: don't
# this is terrible!

function check_sweep_length(sweep)
    @assert length(sweep) == length(sweep.jacobians) == length(sweep.eigvals)
end

function sweeptype(prob::BifurcationProblem,
                   solver::ContinuationSolver,
                   JType::Type = typeof(solver.cache.J),
                   eType::Type = typeof(eigvals_prototype(prob, solver.cache)),
                   )
    tkind = typeof(timekind(solver.cache))
    S = eltype(solver.sol.sweeps)
    uType = eltype(S)
    ptType = point_type_type(prob)
    siType = SpecialPointInterval{tkind, ptType, uType, JType}
    spType = SpecialPoint{tkind, ptType, uType, JType}
    return BifurcationSweep{tkind, ptType, S, JType, eType, siType, spType}
end

function BifurcationSweep(super::ContinuationSweep, solver::ContinuationSolver)
    return sweeptype(solver.prob, solver)(super)
end

struct BifurcationSolution{S <: ContinuationSolution,
                           W <: BifurcationSweep}
    super::S
    sweeps::Vector{W}
end

BifurcationSolution(super::ContinuationSolution,
                    SweesiType::Type{<: BifurcationSweep}) =
    BifurcationSolution(super, SweesiType[])

TimeKind(::Type{<: BifurcationSolution{_, W}}) where {_, W} = TimeKind(W)
