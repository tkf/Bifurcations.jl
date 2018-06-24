module Examples

include("pitchfork.jl")
include("transcritical.jl")
include("calcium.jl")
include("predator_prey.jl")

const PROBLEMS = [
    Pitchfork.prob,
    Transcritical.prob,
    Calcium.prob,
    PredatorPrey.prob,
]
end  # module
