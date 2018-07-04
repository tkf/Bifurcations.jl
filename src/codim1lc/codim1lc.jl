module Codim1LimitCycle  # TODO: merge it to Codim1?  Should I?  Maybe
                         # just rename Codim1 to Codim1FixedPoint?

using ..Continuations: AbstractContinuationProblem, AbstractProblemCache
import ..Continuations: get_prob_cache, get_u0, residual!, residual_jacobian!,
    residual, isindomain

# TimeKind trait:
using ..BifurcationsBase: timekind, Continuous, Discrete
import ..BifurcationsBase: TimeKind

# StateKind trait:
using ..BifurcationsBase: statekind, MutableState, ImmutableState
import ..BifurcationsBase: StateKind

import ..BifurcationsBase: BifurcationProblem

using Setfield: Lens

include("problem.jl")
include("solver.jl")
include("factories.jl")

end  # module
