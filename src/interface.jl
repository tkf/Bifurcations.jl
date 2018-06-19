using .Codim1: Codim1Solver
using .Continuations: ContinuationOptions
import .Continuations: init

init(prob::FixedPointBifurcationProblem; kwargs...) =
    Codim1Solver(prob, ContinuationOptions(; kwargs...))
