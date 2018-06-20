module Bifurcations

# Re-export methods from DifferentialEquations extended here:
export init, solve, solve!, step!
using DiffEqBase: init, solve, solve!, step!
# see: continuations/interface.jl

include("continuations/continuations.jl")
using .Continuations: AbstractContinuationProblem, AbstractProblemCache
import .Continuations: get_prob_cache, get_u0, residual!, residual_jacobian!,
    residual, isindomain
const _C = AbstractProblemCache

include("codim1/codim1.jl")
using .Codim1: timekind, Continuous, Discrete
import .Codim1: TimeKind

include("fixedpoint.jl")
include("interface.jl")
include("diffeq.jl")
include("examples/examples.jl")

# using Requires
# @require RecipesBase include("plotting.jl")
include("plotting.jl")

end # module