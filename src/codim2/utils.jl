"""
    VarDims

A type for calculating conversions between dimensions of the dynamical
system and the augmented system (continuation problem).

| Bifurcation | Variables  | Equations  | `eigvec_dim` | `eigvec_eltype` |
| ----------- | ---------- | ---------- | ------------ | --------------- |
| Saddle-node | ``2n + 2`` | ``2n + 1`` | ``n``        | `<: Real`       |
| Hopf        | ``3n + 3`` | ``3n + 2`` | ``2n``       | `<: Complex`    |

where ``n`` is `ds_dim`, the number of variables is `dim(dom(H))`, and
the number of equations is `dim(codom(H))`.

# Fields
- `ds_dim::Int`: Dimension of the dynamical system.
- `eigvec_dim::Int`: Number of reals required to represent a eigenvector.
- `eigvec_eltype::Type`
"""
struct VarDims
    ds_dim::Int
    eigvec_dim::Int
    eigvec_eltype::Type
end

"""
    augsys_dim(d::VarDims) :: Int

Number of variables `dim(dom(H))` of the augmented system.
"""
augsys_dim(d::VarDims) =
    d.ds_dim + d.eigvec_dim + 2 + (d.eigvec_eltype <: Complex)
eigvec_range(d::VarDims) = d.ds_dim + (1:d.eigvec_dim)
eigval_index(d::VarDims) = d.ds_dim + d.eigvec_dim + 1

function dims_from_augsys(as_dim::Int, T::Type)
    if T <: Real
        ds_dim, r = divrem(as_dim - 2, 2)
        eigvec_dim = ds_dim
    else
        @assert T <: Complex
        ds_dim, r = divrem(as_dim - 3, 3)
        eigvec_dim = 2 * ds_dim
    end
    return _VarDims(ds_dim, eigvec_dim, T, as_dim, r)
end

function _VarDims(ds_dim, eigvec_dim, T, as_dim, r)
    @assert r == 0
    d = VarDims(ds_dim, eigvec_dim, T)
    @assert augsys_dim(d) == as_dim
    return d
end

function dims_from_augsys(as_dom_dim::Int, ::SaddleNodeCont)
    ds_dim, r = divrem(as_dom_dim - 2, 2)
    eigvec_dim = ds_dim
    eigvec_eltype = Real  # TODO: don't
    return _VarDims(ds_dim, eigvec_dim, eigvec_eltype, as_dom_dim, r)
end

function dims_from_augsys(as_dom_dim::Int, ::HopfCont)
    ds_dim, r = divrem(as_dom_dim - 3, 3)
    eigvec_dim = 2 * ds_dim
    eigvec_eltype = Complex  # TODO: don't
    return _VarDims(ds_dim, eigvec_dim, eigvec_eltype, as_dom_dim, r)
end

dims_from_augsys(J::AbstractMatrix, ckind::ContinuationKind) =
    dims_from_augsys(size(J, 2), ckind)

# ds_dim(x::AbstractVector) = ds_dim(length(x))
# ds_dim(var_dim) = var_dim ÷ 2 - 1
eq_dim(N) = 2N + 1
var_dim(N) = 2N + 2

cat_outputs(H1, H2, H3) = vcat(H1, H2, H3)
cat_outputs(H1::SVector, H2::SVector, H3::Number) :: SVector =
    vcat(H1, H2, SVector(H3))

function _as_reals(S, T)
    values = [:($part(v[$i])) for i in 1:S for part in [:real, :imag]]
    quote
        SVector{$(2S), $T}($(values...))
    end
end

as_reals(v::AbstractVector{<: Real}) = v
as_reals(v::Vector{<: Complex{T}}) where T = reinterpret(T, v)

# manually doing this since reinterpret does not work with SVector:
@generated as_reals(v::SVector{S, Complex{T}}) where {S, T} =
    _as_reals(S, T)

as_reals(s::Real) = SVector(s)
as_reals(s::Complex) = SVector(real(s), imag(s))
