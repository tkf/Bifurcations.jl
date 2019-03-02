
<a id='Interface-1'></a>

# Interface

<a id='Bifurcations.BifurcationsBase.BifurcationProblem' href='#Bifurcations.BifurcationsBase.BifurcationProblem'>#</a>
**`Bifurcations.BifurcationsBase.BifurcationProblem`** &mdash; *Type*.



```
BifurcationProblem(point::AbstractSpecialPoint,
                   solver::Codim1Solver,
                   param_axis2::Lens,
                   t2_domain::Tuple)
```

Construct codimension-2 bifurcation problem given a bifurcation `point`.


<a target='_blank' href='https://github.com/tkf/Bifurcations.jl/blob/b17abe3e46f2107dfc63cc3c223adbe86435234e/src/codim2/problem.jl#L9-L17' class='documenter-source'>source</a><br>


```
BifurcationProblem(ode_or_map::AbstractODEProblem,
                   param_axis::Lens,
                   t_domain::Tuple;
                   <keyword arguments>)
```

**Arguments**

  * `ode_or_map`: An `ODEProblem` or `DiscreteProblem`.
  * `param_axis :: Lens`: The lens to set/get a parameter of `ode_or_map`.
  * `t_domain :: Tuple`: A pair of numbers specifying the lower and upper bound for `param_axis`.


<a target='_blank' href='https://github.com/tkf/Bifurcations.jl/blob/b17abe3e46f2107dfc63cc3c223adbe86435234e/src/diffeq.jl#L34-L45' class='documenter-source'>source</a><br>

