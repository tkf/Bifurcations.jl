module TestUtils
using Base.Test
using Compat

using StaticArrays: SVector

using Bifurcations.Codim2: cast_container, as_reals

@testset "cast_container" begin
    # No cast if it's already of the correct type:
    for v in [ones(3),
              SVector(1.0, 2.0),
              SVector(1 + 2im, 3 + 4im)]
        @test cast_container(typeof(v), v) === v
    end

    v = [1.0, 2.0]
    @test cast_container(SVector{2, Float64}, v) === SVector{2, Float64}(v)

    v = [1.0 + 2.0im, 3.0 + 4.0im]
    @test cast_container(SVector{2, Float64}, v) === SVector{2, ComplexF64}(v)
end

@testset "as_reals" begin
    # No convert if it's already a real vector;
    for v in [ones(3),
              SVector(1.0, 2.0)]
        @test as_reals(v) === v
    end

    vc = [1.0 + 2.0im, 3.0 + 4.0im]
    vr = [1.0, 2.0, 3.0, 4.0]
    @test as_reals(vc) == vr
    @test as_reals(SVector{2, ComplexF64}(vc)) == SVector{4, Float64}(vr)
end

end  # module
