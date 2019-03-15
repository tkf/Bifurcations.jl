
<a id='Internals-1'></a>

# Internals


<a id='Continuation-problem-1'></a>

## Continuation problem

<a id='Bifurcations.Continuations.AbstractContinuationProblem' href='#Bifurcations.Continuations.AbstractContinuationProblem'>#</a>
**`Bifurcations.Continuations.AbstractContinuationProblem`** &mdash; *Type*.



Definition of continuation problem.

Numerical continuation algorithms find curves in $\mathbb R^{N}$ implicitly defined by

$$
H(u) = 0
$$

where $H: \mathbb R^{N} \to \mathbb R^{N-1}$.

A continuation problem type (a subtype of `AbstractContinuationProblem`) defines problems for such algorithms to solve by providing how to compute:

  * $H(u)$: [`residual`](internals.md#Bifurcations.Continuations.residual), [`residual!`](internals.md#Bifurcations.Continuations.residual!)
  * its derivative $\partial H /  \partial u$: [`residual_jacobian!`](internals.md#Bifurcations.Continuations.residual_jacobian!)
  * an initial guess $u_0$: [`get_u0`](internals.md#Bifurcations.Continuations.get_u0)
  * and computation cache: [`get_prob_cache`](internals.md#Bifurcations.Continuations.get_prob_cache).


<a target='_blank' href='https://github.com/tkf/Bifurcations.jl/blob/828218782b02b8c931806f1ffdcfa840183df46b/src/continuations/interface.jl#L1-L21' class='documenter-source'>source</a><br>

<a id='Bifurcations.Continuations.AbstractProblemCache' href='#Bifurcations.Continuations.AbstractProblemCache'>#</a>
**`Bifurcations.Continuations.AbstractProblemCache`** &mdash; *Type*.



```
AbstractProblemCache{P <: AbstractContinuationProblem}
```

Cache for computing $H$ and its Jacobian.


<a target='_blank' href='https://github.com/tkf/Bifurcations.jl/blob/828218782b02b8c931806f1ffdcfa840183df46b/src/continuations/interface.jl#L24-L28' class='documenter-source'>source</a><br>

<a id='Bifurcations.Continuations.residual' href='#Bifurcations.Continuations.residual'>#</a>
**`Bifurcations.Continuations.residual`** &mdash; *Function*.



```
residual(u, cache) ↦ H
```

Compute $H(u) \in \mathbb R^{N - 1}$ (aka in-place computation). The definition of $H$ is specified by `cache`.

The name `residual` of the function is came from the problem we are to solve: i.e., find the set of $u$ such that $H(u) = 0$.  Thus, the vector returned by `residual(u, cache)` is considered to be a residual.

**Arguments**

  * `u::AbstractVector` (size: `(N,)`)
  * `cache::AbstractProblemCache`


<a target='_blank' href='https://github.com/tkf/Bifurcations.jl/blob/828218782b02b8c931806f1ffdcfa840183df46b/src/continuations/interface.jl#L31-L45' class='documenter-source'>source</a><br>

<a id='Bifurcations.Continuations.residual!' href='#Bifurcations.Continuations.residual!'>#</a>
**`Bifurcations.Continuations.residual!`** &mdash; *Function*.



```
residual!(H, u, cache) ↦ H
```

Compute $H(u) \in \mathbb R^{N - 1}$ and store the result in `H` (aka out-of-place computation). The definition of $H$ is specified by `cache`.

See also: [`residual`](internals.md#Bifurcations.Continuations.residual)

**Arguments**

  * `H::AbstractVector` (size: `(N - 1,)`)
  * `u::AbstractVector` (size: `(N,)`)
  * `cache::AbstractProblemCache`


<a target='_blank' href='https://github.com/tkf/Bifurcations.jl/blob/828218782b02b8c931806f1ffdcfa840183df46b/src/continuations/interface.jl#L48-L61' class='documenter-source'>source</a><br>

<a id='Bifurcations.Continuations.residual_jacobian!' href='#Bifurcations.Continuations.residual_jacobian!'>#</a>
**`Bifurcations.Continuations.residual_jacobian!`** &mdash; *Function*.



```
residual_jacobian!(H, J, u, cache) ↦ (H, J)
```

Compute $H(u)$ and its Jacobian $\partial H /  \partial u$.

**Arguments**

  * `H::AbstractVector` (size: `(N - 1,)`) $= H(u)$
  * `J::AbstractMatrix` (size: `(N - 1, N)`) $= \partial H /  \partial u$
  * `u::AbstractVector` (size: `(N,)`)
  * `cache::AbstractProblemCache`


<a target='_blank' href='https://github.com/tkf/Bifurcations.jl/blob/828218782b02b8c931806f1ffdcfa840183df46b/src/continuations/interface.jl#L64-L74' class='documenter-source'>source</a><br>

<a id='Bifurcations.Continuations.isindomain' href='#Bifurcations.Continuations.isindomain'>#</a>
**`Bifurcations.Continuations.isindomain`** &mdash; *Function*.



```
isindomain(u, cache) :: Bool
```

**Arguments**

  * `u::AbstractVector` (size: `(N,)`)
  * `cache::AbstractProblemCache`


<a target='_blank' href='https://github.com/tkf/Bifurcations.jl/blob/828218782b02b8c931806f1ffdcfa840183df46b/src/continuations/interface.jl#L77-L83' class='documenter-source'>source</a><br>

<a id='Bifurcations.Continuations.get_prob_cache' href='#Bifurcations.Continuations.get_prob_cache'>#</a>
**`Bifurcations.Continuations.get_prob_cache`** &mdash; *Function*.



```
get_prob_cache(prob::AbstractContinuationProblem) :: AbstractProblemCache
```


<a target='_blank' href='https://github.com/tkf/Bifurcations.jl/blob/828218782b02b8c931806f1ffdcfa840183df46b/src/continuations/interface.jl#L86-L88' class='documenter-source'>source</a><br>

<a id='Bifurcations.Continuations.get_u0' href='#Bifurcations.Continuations.get_u0'>#</a>
**`Bifurcations.Continuations.get_u0`** &mdash; *Function*.



```
get_u0(prob::AbstractContinuationProblem) ↦ u0
```


<a target='_blank' href='https://github.com/tkf/Bifurcations.jl/blob/828218782b02b8c931806f1ffdcfa840183df46b/src/continuations/interface.jl#L93-L95' class='documenter-source'>source</a><br>


<a id='Continuation-algorithm-1'></a>

## Continuation algorithm

<a id='Bifurcations.Continuations.ContinuationCache' href='#Bifurcations.Continuations.ContinuationCache'>#</a>
**`Bifurcations.Continuations.ContinuationCache`** &mdash; *Type*.



Cache for Euler-Newton continuation method.

See [`AbstractContinuationProblem`](internals.md#Bifurcations.Continuations.AbstractContinuationProblem) for the mathematical setup.

**Fields**

  * `prob_cache`
  * `u` (size: `(N,)`)
  * `H` (size: `(N - 1,)`) $= H(u)$
  * `J` (size: `(N - 1, N)`) $= \partial H / \partial u$
  * `Q` (size: `(N - 1, N)`): temporary array for the QR decomposition
  * `h::Real`: step size
  * `direction::Int`: +1 or -1
  * `corrector_success::Bool`
  * `adaptation_success::Bool`
  * `simple_bifurcation::Bool`


<a target='_blank' href='https://github.com/tkf/Bifurcations.jl/blob/828218782b02b8c931806f1ffdcfa840183df46b/src/continuations/euler_newton.jl#L12-L28' class='documenter-source'>source</a><br>


<a id='Bifurcation-problem-1'></a>

## Bifurcation problem

<a id='Bifurcations.FixedPointBifurcationProblem' href='#Bifurcations.FixedPointBifurcationProblem'>#</a>
**`Bifurcations.FixedPointBifurcationProblem`** &mdash; *Type*.



Fixed point bifurcation problem.

See also: [`AbstractContinuationProblem`](internals.md#Bifurcations.Continuations.AbstractContinuationProblem)

**Fields**

  * `homotopy::Function`: A function to compute $H(x, t)$ where $H$ is a homotopy $H: \mathbb R^N \times \mathbb R \to \mathbb R^N$. Function `homotopy` must be callable in one of the following form: `homotopy(x, p, t) ↦ H` (return `H`) for mutable state type or `homotopy(H, x, p, t)` (mutate `H`) for immutable state type.
  * `homotopy_jacobian::Union{Function, Nothing}`: A function to compute $H(x, t)$ and its Jacobian $J = \partial H / \partial (x, t) \in \mathbb R^{N \times (N+1)}$. Function `homotopy_jacobian` must be callable in one of the following form: `homotopy_jacobian(x, p, t) ↦ (H, J)` (return `(H, J)`) or `homotopy_jacobian(H, J, x, p, t)` (mutate `H` and `J`).
  * `u0::Union{AbstractArray, Real}`: Initial state.
  * `t0::Real`: Initial parameter.
  * `t_domain::Tuple{<:Real, <:Real}`: Range of the parameter.
  * `phase_space::Tuple{typeof(u0), typeof(u0)}`: A pair of lower and upper bound of the phase space.  Default is unbounded.
  * `p`: Model parameter (constants).


<a target='_blank' href='https://github.com/tkf/Bifurcations.jl/blob/828218782b02b8c931806f1ffdcfa840183df46b/src/fixedpoint.jl#L8-L32' class='documenter-source'>source</a><br>

