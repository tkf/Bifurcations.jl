tkindstr(x) = tkindstr(timekind(x))
tkindstr(::Discrete) = "Discrete"
tkindstr(::Continuous) = "Continuous"

function Base.show(io::IO, sweep::Codim1Sweep)
    super = as(sweep, ContinuationSweep)
    print(io, "Codim1Sweep <", tkindstr(sweep), ">")

    n_all = length(sweep)
    n_sb = length(super.simple_bifurcation)
    n_sp = length(sweep.special_points)
    if get(io, :compact, false)
        println(io, " ", n_sb + n_sp, "/", n_all, " special/points")
    else
        println(io)
        println(io, "# points             : ", n_all)
        println(io, "# simple bifurcations: ", n_sb)
        println(io, "# special points     : ", n_sp)
    end
end

function show_solution_info(io::IO, sol::Codim1Solution)
    super = as(sol, ContinuationSolution)
    if isempty(sol.sweeps)
        println(io, " no sweeps")
        return
    end
    n_all = sum(length, super.sweeps)
    n_sb = sum(length(s.simple_bifurcation) for s in super.sweeps)
    n_sp = sum(length(s.special_points) for s in sol.sweeps)
    if get(io, :compact, false)
        println(io, " ", n_sb + n_sp, "/", n_all, " special/points")
    else
        println(io)
        println(io, "# sweeps             : ", length(sol.sweeps))
        println(io, "# points             : ", n_all)
        println(io, "# simple bifurcations: ", n_sb)
        println(io, "# special points     : ", n_sp)
    end
end

function Base.show(io::IO, sol::Codim1Solution)
    tkind = isempty(sol.sweeps) ? "?" : tkindstr(sol.sweeps[1]) # FIXME
    print(io, "Codim1Solution <", tkind, ">")
    show_solution_info(io, sol)
end

function Base.show(io::IO, solver::Codim1Solver)
    sol = solver.sol
    tkind = isempty(sol.sweeps) ? "?" : tkindstr(sol.sweeps[1]) # FIXME
    print(io, "Codim1Solver <", tkind, ">")
    show_solution_info(io, sol)
end
