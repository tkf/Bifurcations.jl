module TestUtils
using Base.Test
using Compat

using StaticArrays: SVector

using Bifurcations.Codim2: cast_container, as_reals, _ds_eigvec

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

@testset "_ds_eigvec" begin
    @test _ds_eigvec(Float64, [0, 0, 1, 1, 2, 2]) == [1, 1]
    @test _ds_eigvec(ComplexF64, [0, 0, 1, 1, 1, 1, 3, 2, 2]) ==
        [1 + 1im, 1 + 1im]

    @test _ds_eigvec(Float64, SVector(0, 0, 1, 1, 2, 2)) ==
        SVector(1, 1)
    @test _ds_eigvec(ComplexF64, SVector(0, 0, 1, 1, 1, 1, 3, 2, 2)) ==
        SVector(1 + 1im, 1 + 1im)
end

end  # module
