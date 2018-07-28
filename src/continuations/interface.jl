"""
Definition of continuation problem.

Numerical continuation algorithms find curves in ``\\mathbb R^{N}``
implicitly defined by

```math
H(u) = 0
```

where ``H: \\mathbb R^{N} \\to \\mathbb R^{N-1}``.

A continuation problem type (a subtype of
`AbstractContinuationProblem`) defines problems for such algorithms to
solve by providing how to compute:

* ``H(u)``: [`residual`](@ref), [`residual!`](@ref)
* its derivative ``\\partial H /  \\partial u``: [`residual_jacobian!`](@ref)
* an initial guess ``u_0``: [`get_u0`](@ref)
* and computation cache: [`get_prob_cache`](@ref).
"""
abstract type AbstractContinuationProblem end

"""
    AbstractProblemCache{P <: AbstractContinuationProblem}

Cache for computing ``H`` and its Jacobian.
"""
abstract type AbstractProblemCache{P <: AbstractContinuationProblem} end

"""
    residual(u, cache) ↦ H

Compute ``H(u) \\in \\mathbb R^{N - 1}`` (aka in-place computation).
The definition of ``H`` is specified by `cache`.

The name `residual` of the function is came from the problem we are to
solve: i.e., find the set of ``u`` such that ``H(u) = 0``.  Thus, the
vector returned by `residual(u, cache)` is considered to be a
residual.

# Arguments
- `u::AbstractVector` (size: `(N,)`)
- `cache::AbstractProblemCache`
"""
function residual end  # TODO: remove it (make it internal)

"""
    residual!(H, u, cache) ↦ H

Compute ``H(u) \\in \\mathbb R^{N - 1}`` and store the result in `H`
(aka out-of-place computation).
The definition of ``H`` is specified by `cache`.

See also: [`residual`](@ref)

# Arguments
- `H::AbstractVector` (size: `(N - 1,)`)
- `u::AbstractVector` (size: `(N,)`)
- `cache::AbstractProblemCache`
"""
function residual! end

"""
    residual_jacobian!(H, J, u, cache) ↦ (H, J)

Compute ``H(u)`` and its Jacobian ``\\partial H /  \\partial u``.

# Arguments
- `H::AbstractVector` (size: `(N - 1,)`) ``= H(u)``
- `J::AbstractMatrix` (size: `(N - 1, N)`) ``= \\partial H /  \\partial u``
- `u::AbstractVector` (size: `(N,)`)
- `cache::AbstractProblemCache`
"""
function residual_jacobian! end

"""
    isindomain(u, cache) :: Bool

# Arguments
- `u::AbstractVector` (size: `(N,)`)
- `cache::AbstractProblemCache`
"""
function isindomain end

"""
    get_prob_cache(prob::AbstractContinuationProblem) :: AbstractProblemCache
"""
function get_prob_cache end
# TODO: Rename AbstractProblemCache to ProblemCache and use it instead
# of get_prob_cache.

"""
    get_u0(prob::AbstractContinuationProblem) ↦ u0
"""
function get_u0 end

"""
    after_correction!(prob_cache, u)
"""
function after_correction!(::AbstractProblemCache, ::Any)
end

abstract type AbstractContinuationCache{PC <: AbstractProblemCache} end
