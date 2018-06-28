module BifurcationsBase
using Compat
using ..Continuations: AbstractContinuationCache, AbstractProblemCache,
    AbstractContinuationSolver
include("timekind.jl")
include("statekind.jl")
include("contkind.jl")
include("problem.jl")
include("solution.jl")
include("solver.jl")
include("interface.jl")
include("tools.jl")
include("display.jl")
end  # module
