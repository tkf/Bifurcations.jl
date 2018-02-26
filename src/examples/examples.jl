module Examples

include("pitchfork.jl")
include("transcritical.jl")
include("calcium.jl")

const PROBLEMS = [
    Pitchfork.prob,
    Transcritical.prob,
    Calcium.prob,
]
end  # module
