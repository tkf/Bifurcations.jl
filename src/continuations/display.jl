function print_header(io::IO, x)
    print(io, nameof(typeof(x)))
end

function Base.show(io::IO, sweep::ContinuationSweep)
    print_header(io, sweep)

    n_all = length(sweep)
    n_sb = length(sweep.simple_bifurcation)
    println(io, " ", n_sb, "/", n_all, " bifurcations/points")
end

function show_solution_info(io::IO, sol::ContinuationSolution)
    if isempty(sol.sweeps)
        println(io, " no sweeps")
        return
    end
    n_all = sum(length, sol.sweeps)
    n_sb = sum(length(s.simple_bifurcation) for s in sol.sweeps)
    println(io, " ", n_sb, "/", n_all, " bifurcations/points")
end

function Base.show(io::IO, sol::ContinuationSolution)
    print_header(io, sol)
    show_solution_info(io, sol)
end

function Base.show(io::IO, solver::ContinuationSolver)
    print_header(io, solver)
    show_solution_info(io, solver.sol)
end
