module Continuations
using Compat
import DiffEqBase: solve, solve!, init, step!
include("utils.jl")
include("base.jl")
include("euler_newton.jl")
include("solution.jl")
include("solver.jl")
include("zero_point.jl")
include("branching.jl")
include("display.jl")
include("diagnosis.jl")
end
