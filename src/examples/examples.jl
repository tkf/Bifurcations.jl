module Examples

include("pitchfork.jl")
include("calcium.jl")

const PROBLEMS = [
    Pitchfork.prob,
    Calcium.prob,
]
end  # module
