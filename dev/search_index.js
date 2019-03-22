var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#Bifurcations.jl-documentation-1",
    "page": "Home",
    "title": "Bifurcations.jl documentation",
    "category": "section",
    "text": "Pages = [\n    \"examples.md\",\n    \"api.md\",\n    \"internals.md\",\n]\nDepth = 2"
},

{
    "location": "api/#",
    "page": "API",
    "title": "API",
    "category": "page",
    "text": ""
},

{
    "location": "api/#Bifurcations.BifurcationsBase.BifurcationProblem",
    "page": "API",
    "title": "Bifurcations.BifurcationsBase.BifurcationProblem",
    "category": "type",
    "text": "BifurcationProblem(point::AbstractSpecialPoint,\n                   solver::Codim1Solver,\n                   param_axis2::Lens,\n                   t2_domain::Tuple)\n\nConstruct codimension-2 bifurcation problem given a bifurcation point.\n\n\n\n\n\nBifurcationProblem(ode_or_map::AbstractODEProblem,\n                   param_axis::Lens,\n                   t_domain::Tuple;\n                   <keyword arguments>)\n\nArguments\n\node_or_map: An ODEProblem or DiscreteProblem.\nparam_axis :: Lens: The lens to set/get a parameter of ode_or_map.\nt_domain :: Tuple: A pair of numbers specifying the lower and upper bound for param_axis.\n\n\n\n\n\n"
},

{
    "location": "api/#Interface-1",
    "page": "API",
    "title": "Interface",
    "category": "section",
    "text": "BifurcationProblem"
},

{
    "location": "examples/#",
    "page": "Examples",
    "title": "Examples",
    "category": "page",
    "text": ""
},

{
    "location": "examples/#Examples-1",
    "page": "Examples",
    "title": "Examples",
    "category": "section",
    "text": "Pages = [\n    \"examples/calcium.md\",\n    \"examples/van_der_pol.md\",\n    \"examples/morris_lecar.md\",\n]"
},

{
    "location": "examples/calcium/#",
    "page": "Calcium channel model",
    "title": "Calcium channel model",
    "category": "page",
    "text": ""
},

{
    "location": "examples/calcium/#Calcium-channel-model-1",
    "page": "Calcium channel model",
    "title": "Calcium channel model",
    "category": "section",
    "text": "Calcium channel model taken from PyDSTool example.  See:pydstool/Tutorial_Calcium.py at master · robclewley/pydstool\nTutorial - PyDSTool Wiki\nBifurcation Analysis · DifferentialEquations.jlUse QuickTypes.@qstruct_fp to define model parameter:using QuickTypes: @qstruct_fp\n\n@qstruct_fp CalciumParam(\n    vl = -60,\n    vca = 120,\n    i = -220.0,\n    gl = 2,\n    gca = 4,\n    c = 20,\n    v1 = -1.2,\n    v2 = 18,\n)\n\nnothing # hideDefine the model as in DifferentialEquations.jl:using Parameters: @unpack\n\nfunction f(u, p::CalciumParam, t)\n    @unpack vl, vca, i, gl, gca, c, v1, v2 = p\n    v = u[1]\n    w = u[2]\n    dv = (i + gl * (vl - v) - gca * 0.5 * (1 + tanh((v-v1)/v2)) * (v-vca)) / c\n    dw = v-w\n    return SVector(dv, dw)\nendCreate an ODEProblem:using DiffEqBase: ODEProblem\nusing StaticArrays: SVector\n\nu0 = SVector(-170.0, -170.0)\ntspan = (0.0, 30.0)  # ignored by Bifurcations.jl\np = CalciumParam()\node = ODEProblem(f, u0, tspan, p)Create a bifurcation problem:using Bifurcations: BifurcationProblem\nusing Setfield: @lens\n\nparam_axis = @lens _.i\nprob = BifurcationProblem(ode, param_axis, (-300.0, 100.0))\nnothing # hideSolve it:using DiffEqBase: init, solve!\n\nsolver = init(prob)\nsolve!(solver)\nsol = solver.solPlot it:using Plots\n\nplt = plot(sol)\nsavefig(plt, \"calcium-1.png\"); nothing # hide(Image: )Find the left Saddle-Node bifurcation point:using Bifurcations: special_intervals\n\npoint_list = sort!(special_intervals(solver), by=p->p.u0[end])\npoint = point_list[1]Numerical continuation of the Saddle-Node bifurcation point:sn_prob = BifurcationProblem(\n    point,\n    solver,\n    (@lens _.gca),\n    (0.0, 8.0),\n)\nsn_solver = init(sn_prob)\nsolve!(sn_solver)Plot the phase diagram:plt2 = plot(sn_solver.sol)\nsavefig(plt2, \"calcium-2.png\"); nothing # hide(Image: )"
},

{
    "location": "examples/van_der_pol/#",
    "page": "Continuation of limit cycles of the van der Pol oscillator",
    "title": "Continuation of limit cycles of the van der Pol oscillator",
    "category": "page",
    "text": ""
},

{
    "location": "examples/van_der_pol/#Continuation-of-limit-cycles-of-the-van-der-Pol-oscillator-1",
    "page": "Continuation of limit cycles of the van der Pol oscillator",
    "title": "Continuation of limit cycles of the van der Pol oscillator",
    "category": "section",
    "text": "using Bifurcations\nusing Bifurcations: LimitCycleProblem\nusing Bifurcations.Examples.DuffingVanDerPol\n\nusing Plots\nusing OrdinaryDiffEq: Tsit5, remakeCreate an ODEProblem and solve it:ode = remake(\n    DuffingVanDerPol.ode,\n    p = DuffingVanDerPol.DuffingVanDerPolParam(\n        d = 0.1,\n    ),\n    u0 = [1.0, 1.8],\n    tspan = (0.0, 90),\n)\nsol = solve(ode, Tsit5())\n\nplt_ode = plot(sol, vars=1, tspan=(70, 90))\nsavefig(plt_ode, \"van_der_pol-ode.png\"); nothing # hide(Image: )Let\'s find a point (approximately) on the limit cycle and its period:using Roots: find_zero\nt0 = find_zero((t) -> sol(t)[1] - 1, (80, 83))\nt1 = find_zero((t) -> sol(t)[1] - 1, (t0 + 3, t0 + 7))\nx0 = sol(t0)\n@assert all(isapprox.(x0, sol(t1); rtol=1e-2))\nx0Then a LimitCycleProblem can be constructed from the ode.num_mesh = 50\ndegree = 5\nt_domain = (0.01, 4.0)  # so that it works with this `num_mesh` / `degree`\nprob = LimitCycleProblem(\n    ode, DuffingVanDerPol.param_axis, t_domain,\n    num_mesh, degree;\n    x0 = x0,\n    l0 = t1 - t0,\n    de_args = [Tsit5()],\n)\nnothing  # hideAs the limit cycle is only approximately specified, solver option start_from_nearest_root = true must be passed to start continuation:solver = init(\n    prob;\n    start_from_nearest_root = true,\n    max_branches = 0,\n)\n@time solve!(solver)By default, plot_state_space plots limit cycles colored by its period:plt_lc = plot_state_space(solver)\nsavefig(plt_lc, \"van_der_pol-lc.png\"); nothing # hide(Image: )"
},

{
    "location": "examples/morris_lecar/#",
    "page": "Modified Morris-Lecar model",
    "title": "Modified Morris-Lecar model",
    "category": "page",
    "text": ""
},

{
    "location": "examples/morris_lecar/#Modified-Morris-Lecar-model-1",
    "page": "Modified Morris-Lecar model",
    "title": "Modified Morris-Lecar model",
    "category": "section",
    "text": "Modified Morris-Lecar model from Dhooge, Govaerts, Kuznetsov (2003):Dhooge, Govaerts, Kuznetsov (2003). Numerical Continuation of Fold Bifurcations of Limit Cycles in MATCONTusing Bifurcations\nusing Bifurcations: special_intervals\nusing Bifurcations.Codim1\nusing Bifurcations.Codim2\nusing Bifurcations.Codim2LimitCycle: FoldLimitCycleProblem\nusing Bifurcations.Examples: MorrisLecar\n\nusing Setfield: @lens\nusing PlotsSolve continuation of the equilibrium point:solver = init(\n    MorrisLecar.make_prob();\n    start_from_nearest_root = true,\n    max_branches = 0,\n    nominal_angle_rad = 2π * (5 / 360),\n)\n@time solve!(solver)Plot equilibriums in (u_1 y)-space:plt1 = plot(solver)\nsavefig(plt1, \"morris_lecar-1.png\"); nothing # hide(Image: )"
},

{
    "location": "examples/morris_lecar/#Start-continuation-of-Hopf-bifurcation-1",
    "page": "Modified Morris-Lecar model",
    "title": "Start continuation of Hopf bifurcation",
    "category": "section",
    "text": "hopf_point, = special_intervals(solver, Codim1.PointTypes.hopf)Solve continuation of the Hopf point:codim2_prob = BifurcationProblem(\n    hopf_point,\n    solver,\n    (@lens _.z),\n    (-1.0, 1.0),\n)\nhopf_solver1 = init(\n    codim2_prob;\n    nominal_angle_rad = 0.01,\n)\n@time solve!(hopf_solver1)"
},

{
    "location": "examples/morris_lecar/#Start-continuation-of-fold-bifurcation-of-limit-cycle-at-Bautin-bifurcation-1",
    "page": "Modified Morris-Lecar model",
    "title": "Start continuation of fold bifurcation of limit cycle at Bautin bifurcation",
    "category": "section",
    "text": "bautin_point, = special_intervals(hopf_solver1, Codim2.PointTypes.bautin)Construct a problem for fold bifurcation of the limit cycle starting at bautin_point:flc_prob = FoldLimitCycleProblem(\n    bautin_point,\n    hopf_solver1;\n    period_bound = (0.0, 14.0),  # see below\n    num_mesh = 120,\n    degree = 4,\n)\nflc_solver = init(\n    flc_prob;\n    start_from_nearest_root = true,\n    max_branches = 0,\n    bidirectional_first_sweep = false,\n    nominal_angle_rad = 2π * (5 / 360),\n    max_samples = 500,\n)\n@time solve!(flc_solver)Plot the limit cycles at fold bifurcation boundaries:plt_state_space = plot_state_space(flc_solver)\nsavefig(plt_state_space, \"morris_lecar-state_space.png\"); nothing # hide(Image: )The continuation was configured to stop just before the period is about to diverge.  Note that stopping at larger period requires larger mesh size.plt_periods = plot(flc_solver, (x=:p1, y=:period))\nsavefig(plt_periods, \"morris_lecar-periods.png\"); nothing # hide(Image: )"
},

{
    "location": "examples/morris_lecar/#Start-continuation-of-Saddle-Node-bifurcation-1",
    "page": "Modified Morris-Lecar model",
    "title": "Start continuation of Saddle-Node bifurcation",
    "category": "section",
    "text": "sn_point, = special_intervals(solver, Codim1.PointTypes.saddle_node)Going back to the original continuation of the equilibrium, let\'s start continuation of one of the saddle-node bifurcation:sn_prob = BifurcationProblem(\n    sn_point,\n    solver,\n    (@lens _.z),\n    (-1.0, 1.0),\n)\nsn_solver = init(\n    sn_prob;\n    nominal_angle_rad = 0.01,\n    max_samples = 1000,\n    start_from_nearest_root = true,\n)\n@time solve!(sn_solver)"
},

{
    "location": "examples/morris_lecar/#Switching-to-continuation-of-Hopf-bifurcation-at-Bogdanov-Takens-bifurcation-1",
    "page": "Modified Morris-Lecar model",
    "title": "Switching to continuation of Hopf bifurcation at Bogdanov-Takens bifurcation",
    "category": "section",
    "text": "hopf_prob2 = BifurcationProblem(\n    special_intervals(sn_solver, Codim2.PointTypes.bogdanov_takens)[1],\n    sn_solver,\n)\nhopf_solver2 = init(hopf_prob2)\n@time solve!(hopf_solver2)"
},

{
    "location": "examples/morris_lecar/#Phase-diagram-1",
    "page": "Modified Morris-Lecar model",
    "title": "Phase diagram",
    "category": "section",
    "text": "plt2 = plot()\nfor s in [hopf_solver1, flc_solver, sn_solver, hopf_solver2]\n    plot!(plt2, s)\nend\nplot!(plt2, ylim=(0.03, 0.15), xlim=(-0.05, 0.2))\n\nsavefig(plt2, \"morris_lecar-2.png\"); nothing # hide(Image: )"
},

{
    "location": "internals/#",
    "page": "Internals",
    "title": "Internals",
    "category": "page",
    "text": ""
},

{
    "location": "internals/#Internals-1",
    "page": "Internals",
    "title": "Internals",
    "category": "section",
    "text": ""
},

{
    "location": "internals/#Bifurcations.Continuations.AbstractContinuationProblem",
    "page": "Internals",
    "title": "Bifurcations.Continuations.AbstractContinuationProblem",
    "category": "type",
    "text": "Definition of continuation problem.\n\nNumerical continuation algorithms find curves in mathbb R^N implicitly defined by\n\nH(u) = 0\n\nwhere H mathbb R^N to mathbb R^N-1.\n\nA continuation problem type (a subtype of AbstractContinuationProblem) defines problems for such algorithms to solve by providing how to compute:\n\nH(u): residual, residual!\nits derivative partial H   partial u: residual_jacobian!\nan initial guess u_0: get_u0\nand computation cache: get_prob_cache.\n\n\n\n\n\n"
},

{
    "location": "internals/#Bifurcations.Continuations.AbstractProblemCache",
    "page": "Internals",
    "title": "Bifurcations.Continuations.AbstractProblemCache",
    "category": "type",
    "text": "AbstractProblemCache{P <: AbstractContinuationProblem}\n\nCache for computing H and its Jacobian.\n\n\n\n\n\n"
},

{
    "location": "internals/#Bifurcations.Continuations.residual",
    "page": "Internals",
    "title": "Bifurcations.Continuations.residual",
    "category": "function",
    "text": "residual(u, cache) ↦ H\n\nCompute H(u) in mathbb R^N - 1 (aka in-place computation). The definition of H is specified by cache.\n\nThe name residual of the function is came from the problem we are to solve: i.e., find the set of u such that H(u) = 0.  Thus, the vector returned by residual(u, cache) is considered to be a residual.\n\nArguments\n\nu::AbstractVector (size: (N,))\ncache::AbstractProblemCache\n\n\n\n\n\n"
},

{
    "location": "internals/#Bifurcations.Continuations.residual!",
    "page": "Internals",
    "title": "Bifurcations.Continuations.residual!",
    "category": "function",
    "text": "residual!(H, u, cache) ↦ H\n\nCompute H(u) in mathbb R^N - 1 and store the result in H (aka out-of-place computation). The definition of H is specified by cache.\n\nSee also: residual\n\nArguments\n\nH::AbstractVector (size: (N - 1,))\nu::AbstractVector (size: (N,))\ncache::AbstractProblemCache\n\n\n\n\n\n"
},

{
    "location": "internals/#Bifurcations.Continuations.residual_jacobian!",
    "page": "Internals",
    "title": "Bifurcations.Continuations.residual_jacobian!",
    "category": "function",
    "text": "residual_jacobian!(H, J, u, cache) ↦ (H, J)\n\nCompute H(u) and its Jacobian partial H   partial u.\n\nArguments\n\nH::AbstractVector (size: (N - 1,)) = H(u)\nJ::AbstractMatrix (size: (N - 1, N)) = partial H   partial u\nu::AbstractVector (size: (N,))\ncache::AbstractProblemCache\n\n\n\n\n\n"
},

{
    "location": "internals/#Bifurcations.Continuations.isindomain",
    "page": "Internals",
    "title": "Bifurcations.Continuations.isindomain",
    "category": "function",
    "text": "isindomain(u, cache) :: Bool\n\nArguments\n\nu::AbstractVector (size: (N,))\ncache::AbstractProblemCache\n\n\n\n\n\n"
},

{
    "location": "internals/#Bifurcations.Continuations.get_prob_cache",
    "page": "Internals",
    "title": "Bifurcations.Continuations.get_prob_cache",
    "category": "function",
    "text": "get_prob_cache(prob::AbstractContinuationProblem) :: AbstractProblemCache\n\n\n\n\n\n"
},

{
    "location": "internals/#Bifurcations.Continuations.get_u0",
    "page": "Internals",
    "title": "Bifurcations.Continuations.get_u0",
    "category": "function",
    "text": "get_u0(prob::AbstractContinuationProblem) ↦ u0\n\n\n\n\n\n"
},

{
    "location": "internals/#Continuation-problem-1",
    "page": "Internals",
    "title": "Continuation problem",
    "category": "section",
    "text": "Bifurcations.Continuations.AbstractContinuationProblem\nBifurcations.Continuations.AbstractProblemCache\nBifurcations.Continuations.residual\nBifurcations.Continuations.residual!\nBifurcations.Continuations.residual_jacobian!\nBifurcations.Continuations.isindomain\nBifurcations.Continuations.get_prob_cache\nBifurcations.Continuations.get_u0"
},

{
    "location": "internals/#Bifurcations.Continuations.ContinuationCache",
    "page": "Internals",
    "title": "Bifurcations.Continuations.ContinuationCache",
    "category": "type",
    "text": "Cache for Euler-Newton continuation method.\n\nSee AbstractContinuationProblem for the mathematical setup.\n\nFields\n\nprob_cache\nu (size: (N,))\nH (size: (N - 1,)) = H(u)\nJ (size: (N - 1, N)) = partial H  partial u\nQ (size: (N - 1, N)): temporary array for the QR decomposition\nh::Real: step size\ndirection::Int: +1 or -1\ncorrector_success::Bool\nadaptation_success::Bool\nsimple_bifurcation::Bool\n\n\n\n\n\n"
},

{
    "location": "internals/#Continuation-algorithm-1",
    "page": "Internals",
    "title": "Continuation algorithm",
    "category": "section",
    "text": "Bifurcations.Continuations.ContinuationCache"
},

{
    "location": "internals/#Bifurcations.FixedPointBifurcationProblem",
    "page": "Internals",
    "title": "Bifurcations.FixedPointBifurcationProblem",
    "category": "type",
    "text": "Fixed point bifurcation problem.\n\nSee also: AbstractContinuationProblem\n\nFields\n\nhomotopy::Function: A function to compute H(x t) where H is a homotopy H mathbb R^N times mathbb R to mathbb R^N. Function homotopy must be callable in one of the following form: homotopy(x, p, t) ↦ H (return H) for mutable state type or homotopy(H, x, p, t) (mutate H) for immutable state type.\nhomotopy_jacobian::Union{Function, Nothing}: A function to compute H(x t) and its Jacobian J = partial H  partial (x t) in mathbb R^N times (N+1). Function homotopy_jacobian must be callable in one of the following form: homotopy_jacobian(x, p, t) ↦ (H, J) (return (H, J)) or homotopy_jacobian(H, J, x, p, t) (mutate H and J).\nu0::Union{AbstractArray, Real}: Initial state.\nt0::Real: Initial parameter.\nt_domain::Tuple{<:Real, <:Real}: Range of the parameter.\nphase_space::Tuple{typeof(u0), typeof(u0)}: A pair of lower and upper bound of the phase space.  Default is unbounded.\np: Model parameter (constants).\n\n\n\n\n\n"
},

{
    "location": "internals/#Bifurcation-problem-1",
    "page": "Internals",
    "title": "Bifurcation problem",
    "category": "section",
    "text": "Bifurcations.FixedPointBifurcationProblem"
},

]}
