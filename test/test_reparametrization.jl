module TestReparametrization
include("preamble.jl")

using StaticArrays: SVector

using Bifurcations.Examples.Reparametrization: reparametrize, forward, backward
using Bifurcations.Examples.Bautin

ode0 = reparametrize(Bautin.ode)
ode1 = reparametrize(Bautin.ode; seed=1, extra=SVector(-1.0))
ode2 = reparametrize(Bautin.ode; seed=2, extra=SVector(-1.0, +1.0))

@testset begin
    x1 = [0.1, 0.1]
    y1 = forward(x1, ode0.p)
    x2 = backward(y1, ode0.p)
    y2 = forward(x2, ode0.p)
    @test x1 ≈ x2
    @test y1 ≈ y2

    x1 = [0.1, 0.1, 0.1]
    y1 = forward(x1, ode1.p)
    x2 = backward(y1, ode1.p)
    y2 = forward(x2, ode1.p)
    @test x1 ≈ x2
    @test y1 ≈ y2

    x1 = [0.1, 0.1, 0.1, 0.1]
    y1 = forward(x1, ode2.p)
    x2 = backward(y1, ode2.p)
    y2 = forward(x2, ode2.p)
    @test x1 ≈ x2
    @test y1 ≈ y2
end

@testset "Random test dim=$dim" for (dim, ode) in [(2, ode0),
                                                   (3, ode1),
                                                   (4, ode2)]
    for seed in 1:10
        rng = MersenneTwister(seed)
        x1 = randn(rng, dim)
        y1 = forward(x1, ode.p)
        x2 = backward(y1, ode.p)
        y2 = forward(x2, ode.p)
        @test x1 ≈ x2
        @test y1 ≈ y2
    end
end

end  # module
