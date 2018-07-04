module Codim2

using ..ArrayUtils: cast_container

using ..Continuations: AbstractContinuationProblem, AbstractProblemCache
import ..Continuations: get_prob_cache, get_u0, residual!, residual_jacobian!,
    residual, isindomain

# TimeKind trait:
using ..BifurcationsBase: timekind, Continuous, Discrete
import ..BifurcationsBase: TimeKind

# StateKind trait:
using ..BifurcationsBase: statekind, MutableState, ImmutableState
import ..BifurcationsBase: StateKind

# ContinuationKind trait:
using ..BifurcationsBase: contkind, SaddleNodeCont, HopfCont
using ..BifurcationsBase: ContinuationKind

import ..BifurcationsBase: BifurcationProblem

using Setfield: Lens
using StaticArrays: SVector

include("utils.jl")
include("augmented_systems.jl")
include("problem.jl")
include("diffeq.jl")
include("solver.jl")
include("analysis.jl")
include("switch.jl")
include("api.jl")
include("display.jl")
include("diagnosis.jl")

end  # module
