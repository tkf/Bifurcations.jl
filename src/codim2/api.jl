using ..Continuations: ContinuationOptions
import ..Continuations: init

init(prob::Codim2Problem; kwargs...) =
    BifurcationSolver(prob, ContinuationOptions(; kwargs...))
