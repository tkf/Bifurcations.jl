module Continuations
using LinearAlgebra
using StaticArrays: StaticArray, SMatrix, SVector, SArray
import DiffEqBase: solve, solve!, init, step!

using ..ArrayUtils: _similar, _zeros, isalmostzero, zero_if_nan, _lq!, _det,
    _normalize!, bottomrow, popbottomright

include("interface.jl")
include("euler_newton.jl")
include("solution.jl")
include("solver.jl")
include("zero_point.jl")
include("branching.jl")
include("display.jl")
include("diagnosis.jl")
end
