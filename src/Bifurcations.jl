module Bifurcations

# Re-export methods from DifferentialEquations extended here:
export init, solve, solve!, step!
import DiffEqBase: init, solve, solve!, step!
# see: solver.jl

include("continuations/continuations.jl")

# Continuation algorithm interface:
using .Continuations: ContinuationCache, ContinuationOptions,
    predictor_corrector_step!, new_branches!
# TODO: Reconsider what/how `.Continuations` module "exposes"
# continuation algorithm/method/solver interface.  For example, should
# `predictor_corrector_step!` be renamed to `step!`?

# Continuation problem interface:
using .Continuations: AbstractContinuationProblem, AbstractProblemCache
import .Continuations: get_prob_cache, get_u0, residual!, residual_jacobian!,
    residual, isindomain
const _C = AbstractProblemCache

include("utils.jl")

include("solution.jl")
include("solver.jl")

include("fixedpoint.jl")
include("diffeq.jl")
include("examples/examples.jl")

# using Requires
# @require RecipesBase include("plotting.jl")
include("plotting.jl")

end # module
