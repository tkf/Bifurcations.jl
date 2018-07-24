module FDiffUtils

using ForwardDiff
using ForwardDiff: Dual
using Parameters: @with_kw
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

nondual(x::Dual) = nondual(x.value)
nondual(x::AbstractArray) = nondual.(x)
nondual(x) = x

_all(f, x::AbstractArray) = all(f, x)
_all(f, x::Dual) = _all(f, x.value) && all(_all.(Ref(f), x.partials))
_all(f, x) = f(x)
_any(f, x) = ! _all(!f, x)

all_isfinite(x) = _all(isfinite, x)
any_isnan(x) = _any(isnan, x)

@with_kw struct NotFinite <: Exception
    value
    expr = nothing
end

function Base.showerror(io::IO, err::NotFinite)
    if err.expr != nothing
        show(io, expr)
    else
        print(io, "Value")
    end
    println(io, " contains non-finite value(s):")
    show(IOContext(io, :compact => false), err.value)
    println(io)
end

macro assert_finite(expr)
    quote
        x = $(esc(expr))
        if ! all_isfinite(x)
            error(NotFinite(x, $(QuoteNode(expr))))
        end
        x
    end
end

end  # module
