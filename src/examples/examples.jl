module Examples

# Utilities
include("reparametrization.jl")

# Dynamical systems
include("pitchfork.jl")
include("transcritical.jl")
include("bautin.jl")
include("cusp.jl")
include("calcium.jl")
include("predator_prey.jl")
include("bazykin_85.jl")
include("bazykin_khibnik_81.jl")
include("duffing_van_der_pol.jl")
include("morris_lecar.jl")


example_modules() = [
    Pitchfork,
    Transcritical,
    Bautin,
    Cusp,
    Calcium,
    PredatorPrey,
    Bazykin85,
    BazykinKhibnik81,
    DuffingVanDerPol,
    MorrisLecar,
]

examples() = [nameof(ex) => ex for ex in example_modules()]

end  # module

using .Examples: examples
using .Examples.Reparametrization: reparametrize
