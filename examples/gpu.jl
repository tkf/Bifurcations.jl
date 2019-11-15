using Bifurcations
using LinearAlgebra
using OrdinaryDiffEq
using Setfield
using Parameters: @unpack

HAS_CUDA = try
    using CuArrays
    import CUDAnative
    true
catch err
    @error "CuArrays cannot be imported" exception = (err, catch_backtrace())
    global const CuArray = Nothing
    false
end

function phase_dynamics!(du, u, p, t)
    @unpack J, h, g = p
    tanh = u isa CuArray ? CUDAnative.tanh : Base.tanh
    du .= .- u .+ g .* (J * tanh.(u)) .+ h
    return du
end

function jacobian!(M, u, p, t)
    @unpack J, h, g = p
    cosh = u isa CuArray ? CUDAnative.cosh : Base.cosh
    M .= -I + g .* J .* Diagonal(inv.(cosh.(u)) .^ 2)
    return M
end

N = 10000
ode = ODEProblem(
    ODEFunction(phase_dynamics!; jac = jacobian!),
    randn(N),
    (0.0, 50.0),
    (
        J = randn(N, N) ./ .√(N),
        h = randn(N),
        g = 0.5,
    ),
)

ts = solve(ode, Tsit5())
ode0 = remake(ode, u0 = ts.u[end])

param_axis = @lens _.g
prob = BifurcationProblem(ode0, param_axis, (0.5, 1.5))
# prob = BifurcationProblem(ode, param_axis, (0.5, 1.5))

using Bifurcations.Continuations: ContinuationSolver, ContinuationOptions
opts = ContinuationOptions(
    atol = 0.05,
    rtol = 0.05,
    max_branches = 0,
    max_samples = 1,
    nominal_distance = 0.01N,
    # bidirectional_first_sweep = false,
    # start_from_nearest_root = true,
    verbose = true,
)
solver = ContinuationSolver(prob, opts)
# solve(prob; start_from_nearest_root = true)
# sol = @time solve!(solver)
# sol = @time solve!(solver)

if HAS_CUDA
    ode_gpu = remake(
        ode0;
        u0 = cu(ode0.u0),
        p = (
            J = cu(ode0.p.J),
            h = cu(ode0.p.h),
            g = ode0.p.g,
        ),
    )
    prob_gpu = BifurcationProblem(ode_gpu, param_axis, (0.5, 1.5))
    solver_gpu = ContinuationSolver(prob_gpu, opts)
    # solve(prob_gpu; start_from_nearest_root = true)
    # sol_gpu = @time solve!(solver_gpu)
    # sol_gpu = @time solve!(solver_gpu)
end
