module Examples

include("pitchfork.jl")
include("transcritical.jl")
include("bautin.jl")
include("cusp.jl")
include("calcium.jl")
include("predator_prey.jl")
include("bazykin_85.jl")
include("duffing_van_der_pol.jl")

using Compat

example_modules() = [
    Pitchfork,
    Transcritical,
    Bautin,
    Cusp,
    Calcium,
    PredatorPrey,
    Bazykin85,
    DuffingVanDerPol,
]

examples() = [nameof(ex) => ex for ex in example_modules()]

end  # module

using .Examples: examples
