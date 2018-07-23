module TestResidualJacobian
include("preamble.jl")
include("sym_diff.jl")

using StaticArrays: SVector

using .SymDiff: test_residual_jacobian
using Bifurcations.BifurcationsBase: SaddleNodeCont, HopfCont
using Bifurcations.Codim2: NormalizingAS, BackReferencingAS
using Bifurcations.Examples: examples

EXAMPLES = [
    (name, ex) for (name, ex) in examples()
    if try
        ex.make_codim2_prob
        true
    catch
        false
    end
]

@testset "residual_jacobian $name" for (name, ex) in EXAMPLES
    @testset "$(nameof(ck)) $(nameof(ASType)) $(utype)" for
        utype in [:SVector, :Vector],
        ASType in [BackReferencingAS, NormalizingAS],
        ck in [SaddleNodeCont, HopfCont]

        if utype == :SVector
            u0 = SVector(ex.u0...)
        elseif utype == :Vector
            u0 = Vector(ex.u0)
        else
            error("utype: $utype")
        end

        v0 = copy(u0)
        w0 = 0
        @assert ck in (SaddleNodeCont, HopfCont)
        if ck == HopfCont
            v0 = v0 .+ 1im
            w0 = 1.0
        end

        prob = ex.make_codim2_prob(
            u0 = u0,
            v0 = v0,
            w0 = w0,
            augmented_system = ASType(),
        )
        test_residual_jacobian(prob)
    end
end

end  # module
