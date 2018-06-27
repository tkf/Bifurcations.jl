module ArrayUtils

using StaticArrays: SMatrix

_eigvals(A) = eigvals(A)
_eigvals(A::SMatrix) = eigvals(Array(A))
# Workaround: Only hermitian matrices are diagonalizable by *StaticArrays*.

_eig(A) = eig(A)
_eig(A::SMatrix) = eig(Array(A))

end  # module
