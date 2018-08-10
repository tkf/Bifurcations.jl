using Printf: @sprintf

set_if_not(io, key, val) = haskey(io, key) ? io : IOContext(io, key => val)

function print_header(io::IO, x)
    print(io, nameof(typeof(x)))
end

function print_header(io::IO, point::SimpleBifurcationInterval)
    i_sweep = try
        point.sweep.value.i
    catch
        '?'
    end
    h = @sprintf "%.4g" norm(point.h)
    print(io, nameof(typeof(point)),
          " u0=sweeps[$i_sweep].u[$(point.i)]",
          " h=", h,
          " dir=", point.direction)
end

function Base.show(io::IO, point::SimpleBifurcationInterval)
    print_header(io, point)
    println(io)
    if ! get(io, :compact, false)
        io = set_if_not(io, :compact, true)  # reduce number of digits shown
        println(io, "happened between:")
        println(io, "  u0 = ", point.u0)
        println(io, "  u1 = ", point.u1)
    end
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
