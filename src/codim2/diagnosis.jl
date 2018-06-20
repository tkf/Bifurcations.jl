import ..Continuations: find_errors

function find_errors(solver::Codim2Solver)
    super = as(solver, ContinuationSolver)
    errors = find_errors(super)
    # TODO: do some Codim2Solver-specific checks
    return errors
end
