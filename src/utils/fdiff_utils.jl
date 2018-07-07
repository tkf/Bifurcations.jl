module FDiffUtils

using ForwardDiff
using StaticArrays: SVector

function vderiv2(f, x::Number = 0.0)
    return ForwardDiff.jacobian(
        y -> ForwardDiff.jacobian(z -> f(z[1]), y)[:, 1],
        SVector(x),
    )[:, 1]
end

function deriv2(f, x::Number = 0.0)
    return ForwardDiff.derivative(
        a -> ForwardDiff.derivative(f, a),
        x)
end

function deriv3(f, x::Number = 0.0)
    return ForwardDiff.derivative(
        a -> ForwardDiff.derivative(
            b -> ForwardDiff.derivative(f, b),
            a),
        x)
end

end  # module
