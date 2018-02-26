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
    return (Q, UpperTriangular(R))
end

function lq!(Q, A)
    F = lqfact!(A)
    A_mul_B!(F[:Q], eye!(Q))  # L = Matrix(F[:L])[...]; but lesser allocation
    return (LowerTriangular(F[:L]), Q)
end

function lq!(_, A::SMatrix)
    Q, R = qr(A')
    return (LowerTriangular(R'), Q)
end

_A_ldiv_B!(A, B) = A_ldiv_B!(A, B)
_A_ldiv_B!(A, B::SVector) = A \ B
