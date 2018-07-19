module TestBifurcations
using Base.Test

@testset "$file" for file in [
        "test_utils.jl",
        "test_normal_form.jl",
        "test_smoke.jl",
        "test_jacobian.jl",
        "test_residual_jacobian.jl",
        "test_calcium.jl",
        "test_predator_prey.jl",
        "test_bazykin_85.jl",
        "test_bautin.jl",
        "test_reparametrization.jl",
        "test_reparametrized_bautin.jl",
        "test_examples.jl",
        "test_vs_svector.jl",
        "test_duffing_van_der_pol.jl",
        "test_hopf_to_lc.jl",
        ]
    @time include(file)
end

end  # module
