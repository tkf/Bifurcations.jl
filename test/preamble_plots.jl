include("preamble.jl")

import GR  # Workaround Plots.jl world age problem
using Plots

using Bifurcations: plot  # TODO: stop doing this


function nullshow(plt::Plots.Plot)
    nullshow(MIME("image/png"), plt)
end
