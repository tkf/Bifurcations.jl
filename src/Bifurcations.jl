module Bifurcations

# Re-export methods from DifferentialEquations extended here:
export init, solve, solve!, step!
using DiffEqBase: init, solve, solve!, step!
# see: continuations/solver.jl

include("utils/utils.jl")

include("continuations/continuations.jl")
using .Continuations: AbstractContinuationProblem, AbstractProblemCache
import .Continuations: get_prob_cache, get_u0, residual!, residual_jacobian!,
    residual, isindomain

include("base/base.jl")
# TimeKind trait:
using .BifurcationsBase: timekind, Continuous, Discrete
import .BifurcationsBase: TimeKind
# StateKind trait:
using .BifurcationsBase: statekind, MutableState, ImmutableState
import .BifurcationsBase: StateKind

using .BifurcationsBase: BifurcationProblem
export BifurcationProblem

include("codim1/codim1.jl")

include("fixedpoint.jl")
include("api.jl")
include("diffeq.jl")

include("codim1lc/codim1lc.jl")
using .Codim1LimitCycle: LimitCycleProblem

include("codim2/codim2.jl")
include("codim2lc/codim2lc.jl")

include("examples/examples.jl")

# using Requires
# @require RecipesBase include("plotting.jl")
include("plotting.jl")
export plot_state_space, plot_state_space!

end # module
