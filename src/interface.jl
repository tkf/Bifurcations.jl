using .Codim1: Codim1Solver, AbstractSpecialPoint
using .Continuations: ContinuationOptions
import .Continuations: init

init(prob::FixedPointBifurcationProblem; kwargs...) =
    Codim1Solver(prob, ContinuationOptions(; kwargs...))

init(point::AbstractSpecialPoint, args...; kwargs...) =
    init(BifurcationProblem(point, args...); kwargs...)
