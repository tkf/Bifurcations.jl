module TestBifurcations
using Base.Test

@testset "$file" for file in [
        "test_normal_form.jl",
        "test_smoke.jl",
        "test_codim2.jl",
        ]
    @time include(file)
end

end  # module
