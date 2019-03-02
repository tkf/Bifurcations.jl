module PolyUtils

using ForwardDiff

# taken from Jacobi.lagrange
function lagrange(i, x, z)
    nz = length(z)

    l = one(z[1])

    for k = 1:(i-1)
        l = l * (x-z[k]) / (z[i]-z[k])
    end

    for k = (i+1):nz
        l = l * (x-z[k]) / (z[i]-z[k])
    end

    return l
end

dlagrange(i, x, z) = ForwardDiff.derivative(x -> lagrange(i, x, z), x)

end  # module
