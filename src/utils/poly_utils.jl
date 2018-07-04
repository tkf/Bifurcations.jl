module PolyUtils

using Jacobi: lagrange
using ForwardDiff

dlagrange(i, x, z) = ForwardDiff.derivative(x -> lagrange(i, x, z), x)

end  # module
