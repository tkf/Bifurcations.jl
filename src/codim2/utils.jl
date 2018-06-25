using Compat.TypeUtils: typename

struct VarDims
    ds_dim::Int
    eigvec_dim::Int
    eigvec_eltype::Type
end

augsys_dim(d::VarDims) =
    d.ds_dim + d.eigvec_dim + 2 + (d.eigvec_eltype <: Complex)
eigvec_range(d::VarDims) = d.ds_dim + (1:d.eigvec_dim)

function dims_from_augsys(as_dim::Int, T::Type)
    if T <: Real
        ds_dim, r = divrem(as_dim - 2, 2)
        eigvec_dim = ds_dim
    else
        @assert T <: Complex
        ds_dim, r = divrem(as_dim - 3, 3)
        @assert r == 0
        eigvec_dim = 2 * ds_dim
    end
    @assert r == 0
    d = VarDims(ds_dim, eigvec_dim, T)
    @assert augsys_dim(d) == as_dim
    return d
end

# ds_dim(x::AbstractVector) = ds_dim(length(x))
# ds_dim(var_dim) = var_dim ÷ 2 - 1
eq_dim(N) = 2N + 1
var_dim(N) = 2N + 2

function ds_state(x::AbstractArray)
    N = length(x) ÷ 2 - 1
    return @view x[1:N]
end

function sarray_slice_epxr(C, x, indices)
    values = [:($x[$i]) for i in indices]
    return Expr(:call, C, values...)
end

@generated function ds_state(x::SVector{S, T}) where {S, T}
    N = S ÷ 2 - 1
    @assert var_dim(N) == S
    C = :(SVector{$N, T})
    return sarray_slice_epxr(C, :x, 1:N)
end

function output_vars(H)
    N = length(H) ÷ 2
    H1 = @view H[1:N]
    H2 = @view H[N+1:2N]
    H3 = @view H[end]
    return H1, H2, H3
end

cat_outputs(H1, H2, H3) = vcat(H1, H2, H3)
cat_outputs(H1::SVector, H2::SVector, H3::Number) :: SVector =
    vcat(H1, H2, SVector(H3))

cast_container(::Type{A}, v::A) where {A <: AbstractArray} = v

function cast_container(::Type{A}, v) where {A <: AbstractArray}
    constructor = typename(A).wrapper  # e.g., Array
    return constructor(v)
end

@generated function cast_container(::Type{<: SVector{S, <:Number}},
                                   v::AbstractArray{T},
                                   ) where {S, T}
    quote
        SVector{$S, $T}(v)
    end
end

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
