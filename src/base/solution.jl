using ..Continuations: as,
    ContinuationSweep, ContinuationSolution,
    ContinuationCache, ContinuationSolver

abstract type AbstractSpecialPoint{tkind <: TimeKind} end
TimeKind(::Type{<: AbstractSpecialPoint{tkind}}) where tkind = tkind()

struct SpecialPoint{tkind, ptType, uType, JType,
                    } <: AbstractSpecialPoint{tkind}
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
                            } <: AbstractSpecialPoint{tkind}
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
        pType <: SpecialPointInterval{tkind, ptType},
        }
    timekind::tkind
    super::S
    jacobians::Vector{JType}
    eigvals::Vector{eType}
    special_points::Vector{pType}
end
# See: [[../continuations/solution.jl::ContinuationSweep]]

BifurcationSweep{tkind, ptType, S, JType, eType, pType}(super::S) where{
    tkind <: TimeKind,
    ptType,
    S <: ContinuationSweep,
    JType <: AbstractArray,
    eType <: AbstractArray,
    pType <: SpecialPointInterval{tkind, ptType},
} = BifurcationSweep{tkind, ptType, S, JType, eType, pType}(
    tkind(),
    super,
    JType[],
    eType[],
    pType[],
)

TimeKind(::Type{<: BifurcationSweep{tkind}}) where tkind = tkind()

Base.length(sweep::BifurcationSweep) = length(as(sweep, ContinuationSweep))

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
    pType = SpecialPointInterval{tkind, ptType, uType, JType}
    return BifurcationSweep{tkind, ptType, S, JType, eType, pType}
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
                    SweepType::Type{<: BifurcationSweep}) =
    BifurcationSolution(super, SweepType[])

TimeKind(::Type{<: BifurcationSolution{_, W}}) where {_, W} = TimeKind(W)
