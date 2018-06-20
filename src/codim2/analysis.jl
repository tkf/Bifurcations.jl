using ..Codim1: ds_eigvals

function guess_point_type(::Discrete, cache, opts)
    return PointTypes.none
end

function guess_point_type(::Continuous, cache, opts)
    return PointTypes.none
end
