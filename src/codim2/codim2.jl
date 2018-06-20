module Codim2

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
using StaticArrays: SVector

include("base.jl")
include("problem.jl")
include("diffeq.jl")

end  # module
