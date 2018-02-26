module Continuations
import DiffEqBase: solve, solve!, init, step!
include("utils.jl")
include("solution.jl")
include("interface.jl")
include("euler_newton.jl")
include("zero_point.jl")
include("branching.jl")
end
