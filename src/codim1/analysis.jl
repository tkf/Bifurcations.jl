using ..ArrayUtils: _eigvals

function ds_eigvals(::Discrete, J)
    ev = _eigvals(J)
    sort!(ev, by=abs, rev=true)
    return ev
end
# TODO: improve it for StaticArrays

function guess_point_type(::Discrete, cache, eigvals, opts)
    old = abs(cache.eigvals[1]) - 1
    new = abs(eigvals[1]) - 1
    if old * new < 0
        if eigvals[1] > 0
            # eigvals[1] ≈ +1
            return PointTypes.saddle_node
        else
            # eigvals[1] ≈ -1
            return PointTypes.period_doubling
        end
    end
    return PointTypes.none
end

isstable(::Discrete, eigvals) = abs(eigvals[1]) < 1

function ds_eigvals(::Continuous, J)
    ev = _eigvals(J)
    sort!(ev, by=real, rev=true)
    return ev
end

function guess_point_type(::Continuous, cache, eigvals, opts)
    old = real(cache.eigvals[1])
    new = real(eigvals[1])
    if old * new < 0
        e1 = mean((cache.eigvals[1], eigvals[1]))
        if imag(e1) > opts.atol
            return PointTypes.hopf
        else
            return PointTypes.saddle_node
        end
    end
    return PointTypes.none
end

isstable(::Continuous, eigvals) = real(eigvals[1]) < 0
