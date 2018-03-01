using ...Bifurcations: FixedPointBifurcationProblem
import ..Continuations: init

init(prob::FixedPointBifurcationProblem; kwargs...) =
    Codim1Solver(prob, ContinuationOptions(; kwargs...))
