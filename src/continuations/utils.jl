using StaticArrays: SVector, SMatrix, StaticArray, Size, similar_type

_similar(x::AbstractArray, dims...) = similar(x, dims)
_similar(x::StaticArray, dims...) = _zeros(x, dims...)
_zeros(x::AbstractArray{T}, dims...) where T = zeros(x, T, dims)
_zeros(x::StaticArray, dims...) = zeros(similar_type(x, Size(dims)))

function isalmostzero(x, rtol, atol)
    n = norm(x)
    return n <= atol + rtol * n
end

zero_if_nan(x) = isnan(x) ? zero(x) : x

@inline function eye!(A)
    A .= 0
    for i in 1:min(size(A)...)
        @inbounds A[i, i] = 1
    end
    A
end

function qr!(Q, A)
    F = qrfact!(A)
    A_mul_B!(F[:Q], eye!(Q))  # Q = Matrix(F[:Q])[...]; but lesser allocation
    return (Q, UpperTriangular(F[:R]))
end

function qr!(_, A::SMatrix)
    Q, R = qr(A)
    # return (Q, R)
    return (Q, UpperTriangular(R))
end

function lq!(Q, A)
    F = lqfact!(A)
    A_mul_B!(F[:Q], eye!(Q))  # L = Matrix(F[:L])[...]; but lesser allocation
    return (LowerTriangular(F[:L]), Q)
end

function lq!(_, A::SMatrix)
    Q, R = qr(A')
    # return (R', Q)
    return (LowerTriangular(R'), Q)
end

function _normalize!(x)
    normalize!(x)
    return x
end

_normalize!(x::SVector) = x ./ norm(x)

# TODO: use \ instead of A_ldiv_B!
_A_ldiv_B!(A, B::T) where T = _A_ldiv_B!(T, A, B)
_A_ldiv_B!(Y, A, B::T) where T = _A_ldiv_B!(T, Y, A, B)

_A_ldiv_B!(::Type{<: SubArray{T, N, P}}, args...) where {T, N, P} =
    _A_ldiv_B!(P, args...)
_A_ldiv_B!(::Type{<: LowerTriangular{T, P}}, args...) where {T, P} =
    _A_ldiv_B!(P, args...)
_A_ldiv_B!(::Type{<: UpperTriangular{T, P}}, args...) where {T, P} =
    _A_ldiv_B!(P, args...)

_A_ldiv_B!(::Type{<:AbstractArray}, A, B) = A_ldiv_B!(A, B)
_A_ldiv_B!(::Type{<:StaticArray}, A, B) = A \ B
_A_ldiv_B!(::Type{<:Vector}, A::SubArray, B) =  A \ B # TODO: make it in-place
_A_ldiv_B!(::Type{<:AbstractArray}, Y, A, B) = A_ldiv_B!(Y, A, B)
_A_ldiv_B!(::Type{<:StaticArray}, _, A, B) = A \ B
