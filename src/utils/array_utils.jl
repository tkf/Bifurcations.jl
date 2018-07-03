module ArrayUtils

using StaticArrays: SMatrix, SVector

container_array_of(::SVector{S}) where {S} = SVector{S}
container_array_of(::SMatrix{S1, S2}) where {S1, S2} = SMatrix{S1, S2}
container_array_of(::Array{T, N}) where {T, N} = Array{<:Any, N}

_eigvals(A) = eigvals(A)
_eigvals(A::SMatrix) = eigvals(Array(A))
# Workaround: Only hermitian matrices are diagonalizable by *StaticArrays*.

_eig(A) = eig(A)
_eig(A::SMatrix) = eig(Array(A))

end  # module
