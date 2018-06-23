aliasof(x) = nameof(typeof(x))

tkindstr(x) = tkindstr(timekind(x))
tkindstr(::Discrete) = "Discrete"
tkindstr(::Continuous) = "Continuous"

function print_header(io::IO, point::Union{SpecialPoint,
                                           SpecialPointInterval})
    print(io, aliasof(point), " <",
          tkindstr(point), " ",
          point.point_type,
          ">")
end

function print_header(io::IO, ctx::Union{BifurcationSweep,
                                         BifurcationSolution,
                                         BifurcationSolver})
    print(io, aliasof(ctx), " <", tkindstr(ctx), ">")
end

set_if_not(io, key, val) = haskey(io, key) ? io : IOContext(io, key => val)

function Base.show(io::IO, point::SpecialPoint)
    print_header(io, point)
    println(io)
    if ! get(io, :compact, false)
        io = set_if_not(io, :compact, true)  # reduce number of digits shown
        println(io, "u = ", point.u)
    end
end

function Base.show(io::IO, point::SpecialPointInterval)
    print_header(io, point)
    println(io)
    if ! get(io, :compact, false)
        io = set_if_not(io, :compact, true)  # reduce number of digits shown
        println(io, "happened between:")
        println(io, "  u0 = ", point.u0)
        println(io, "  u1 = ", point.u1)
    end
end

function Base.show(io::IO, sweep::BifurcationSweep)
    super = as(sweep, ContinuationSweep)
    print_header(io, sweep)

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

function show_solution_info(io::IO, sol::BifurcationSolution)
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

function Base.show(io::IO, sol::BifurcationSolution)
    print_header(io, sol)
    show_solution_info(io, sol)
end

function Base.show(io::IO, solver::BifurcationSolver)
    print_header(io, solver)
    show_solution_info(io, solver.sol)
end
