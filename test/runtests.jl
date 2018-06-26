module TestBifurcations
using Base.Test

@testset "$file" for file in [
        "test_utils.jl",
        "test_normal_form.jl",
        "test_smoke.jl",
        "test_jacobian.jl",
        "test_calcium.jl",
        "test_predator_prey.jl",
        "test_examples.jl",
        ]
    @time include(file)
end

end  # module
