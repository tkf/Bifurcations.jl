using StaticArrays: SVector, Size, similar_type

_similar(x::AbstractArray, dims...) = similar(x, dims)
_similar(x::SVector, dims...) = zeros(similar_type(x, Size(dims)))

function isalmostzero(x, rtol, atol)
    n = norm(x)
    return n <= atol + rtol * n
end

zero_if_nan(x) = isnan(x) ? zero(x) : x
