module TestBifurcations
using Base.Test

@testset "$file" for file in [
        "test_smoke.jl",
        ]
    @time include(file)
end

end  # module
