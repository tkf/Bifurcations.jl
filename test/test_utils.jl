module TestUtils
using Base.Test
using Compat

using StaticArrays: SVector, SMatrix

using Bifurcations.ArrayUtils: container_array_of, _lq!, canonicalize, nan_
using Bifurcations.Codim2: cast_container, as_reals, _ds_eigvec

@testset "container_array_of" begin
    @test container_array_of([1, 2, 3]) == Vector
    @test container_array_of([1 2; 3 4]) == Matrix
    @test container_array_of(SVector(1.0, 2.0, 3.0)) === SVector{3}
    @test container_array_of(SMatrix{2, 2}(1, 2, 3, 4)) === SMatrix{2, 2}
end

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

@testset "_lq!" begin
    rng = MersenneTwister(0)
    for _ in 1:10, array_type in [Array, SMatrix{3, 3}]
        A = array_type(randn(rng, (3, 3)))

        L0, Q0 = lq(Array(A))

        Q1 = copy(Q0)
        L1, Q1 = _lq!(Q1, copy(A))

        @test L0 ≈ L1
        @test Q0 ≈ Q1
    end
end

@testset "canonicalize" begin
    x_list = Any[
        [1.0, 2],
        [1.0 + 2im, 3 + 4im],
        [1.0 + 0im, 0.0 + 1im],
    ]
    for x in [x1 for x0 in x_list for x1 in [x0, SVector(x0...)]]
        y = canonicalize(x)
        @test typeof(y) === typeof(x)
        @test norm(y) ≈ 1
        @test abs(real(y) ⋅ imag(y)) < 1e-7

        # Test idempotence
        z = canonicalize(y)
        @test typeof(z) === typeof(x)
        @test norm(z) ≈ 1
        @test abs(real(z) ⋅ imag(z)) < 1e-7
    end
end

nanmean(args...) = nan_(mean, args...)
ceil_to_int(x) = ceil(Int, x)

@testset "nan_(mean)" begin
    @testset "$v (-> finite)" for v in [
            10:30,
            collect(-3:-1) .+ 0.0,
            ]
        @test nanmean(v) == mean(v)
        @testset "$f" for f in [
                exp,
                sin,
                ceil_to_int,
                ]
            @test nanmean(f, v) == mean(f, v)
        end
    end
    @testset "$v (-> NaN)" for v in [
            Float64[],
            Int[],
            ]
        @test isnan(nanmean(v))
        # Note: `nanmean(f, v)` throws but that's OK; `mean(f, v)` does so too
    end
    for (x, y) in [
            ([1, NaN, 2, 3], 1:3),
            ([1, NaN, 2, 3], collect(1:3) .+ 0.0),
            ]
        @test nanmean(x) == mean(y)
    end

    A = [
        1   2 3   NaN
        NaN 2 NaN 6
        7   8 9 10
    ]
    @test nanmean(A, 1) == [4  4  6  8]
    @test nanmean(A, 2) == [2 4 8.5]'
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
