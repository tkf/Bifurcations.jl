import ..Continuations: find_errors

function find_errors(solver::Codim1Solver)
    super = as(solver, ContinuationSolver)
    errors = find_errors(super)
    # TODO: do some Codim1Solver-specific checks
    return errors
end
