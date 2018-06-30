module Examples

include("pitchfork.jl")
include("transcritical.jl")
include("calcium.jl")
include("predator_prey.jl")
include("bazykin_85.jl")

const PROBLEMS = [
    Pitchfork.prob,
    Transcritical.prob,
    Calcium.prob,
    PredatorPrey.prob,
    Bazykin85.prob,
]
end  # module
