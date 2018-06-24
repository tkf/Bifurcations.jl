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


function preferred_augsys(point)
    if point.point_type == Codim1.PointTypes.saddle_node
        return NormalizingAS()
    end
    return BackReferencingAS()
end
