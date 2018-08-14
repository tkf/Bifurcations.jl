include("preamble.jl")

import GR  # Workaround Plots.jl world age problem
using Plots



function nullshow(plt::Plots.Plot)
    nullshow(MIME("image/png"), plt)
end
