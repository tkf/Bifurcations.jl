using Compat

using Bifurcations: BifurcationProblem, special_points

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
