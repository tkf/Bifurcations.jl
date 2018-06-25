using ..Continuations: as, ContinuationCache

"""
`AugmentedSystem` subtypes specify how to construct an augmented system.
"""
abstract type AugmentedSystem end

augsys(::T) where T = AugmentedSystem(T)

""" Normalize the eigenvector.  Works only for saddle-node. """
struct NormalizingAS <: AugmentedSystem end

""" Project onto the previous eigenvector (*reference vector*). """
struct BackReferencingAS <: AugmentedSystem end


struct NormalizingASCache end

mutable struct BackReferencingASCache{T}
    v::T
end

get_augsys_cache(prob) = _get_augsys_cache(augsys(prob), prob)
_get_augsys_cache(::NormalizingAS, ::Any) = NormalizingASCache()

set_augsys_cache!(::NormalizingASCache, ::ContinuationCache) = nothing

function set_augsys_cache!(augsys_cache::BackReferencingASCache,
                           cont_cache::ContinuationCache)
    prob = cont_cache.prob_cache.prob  # TODO: interface
    augsys_cache.v = ds_eigvec(prob, cont_cache.u)
end

function set_augsys_cache!(wrapper)
    cont_cache = as(wrapper, ContinuationCache)
    augsys_cache = cont_cache.prob_cache.augsys_cache  # TODO: interface
    set_augsys_cache!(augsys_cache, cont_cache)
end


function preferred_augsys(point)
    if point.point_type == Codim1.PointTypes.saddle_node
        return NormalizingAS()
    end
    return BackReferencingAS()
end
