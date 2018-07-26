using Base.Test
using Compat
using Compat: @info
using Setfield: @lens

using Bifurcations
using Bifurcations: BifurcationProblem, special_points
using Bifurcations.Codim2: cast_container

macro test_nothrow(ex)
    quote
        @test begin
            $(esc(ex))
            true
        end
    end
end

function nullshow(x)
    show(devnull, x)
end

function nullshow(mime::MIME, x)
    show(devnull, mime, x)
end

function smoke_test_solver_show(solver)
    @test_nothrow nullshow(solver)
    @test_nothrow nullshow(solver.sol)

    for sweep in solver.sol.sweeps
        @test_nothrow nullshow(sweep)
        for point in special_points(sweep)
            @test_nothrow nullshow(point)
        end
    end
end

function smoke_test_solver_plot(solver)
    for p in special_points(solver)
        @test_nothrow nullshow(plot(p))
    end
    for sweep in solver.sol.sweeps
        @test_nothrow nullshow(plot(sweep))
    end
    @test_nothrow nullshow(plot(solver.sol))
    @test_nothrow nullshow(plot(solver))
end

_generalize_f(f, A) = function(u::U, p, t) where {T, U <: AbstractArray{T}}
    v = cast_container(A, u)
    H = f(v, p, t)
    return cast_container(U, H)
end

generalized_f(mod) = _generalize_f(mod.f, typeof(mod.u0))


function uniquify_points(point_list; digits=3)
    uniquified = []
    seen = []
    for point in point_list
        y = round.(point.u, digits)
        for x in seen
            if x == y
                @goto seen
            end
        end
        push!(uniquified, point)
        push!(seen, y)

        @label seen
    end
    return uniquified
end
