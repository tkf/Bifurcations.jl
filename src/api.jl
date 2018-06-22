using .BifurcationsBase: BifurcationProblem, BifurcationSolver
using .Codim1: AbstractSpecialPoint
using .Continuations: ContinuationOptions
import .Continuations: init

init(prob::BifurcationProblem; kwargs...) =
    BifurcationSolver(prob, ContinuationOptions(; kwargs...))

init(point::AbstractSpecialPoint, args...; kwargs...) =
    init(BifurcationProblem(point, args...); kwargs...)
