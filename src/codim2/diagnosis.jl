import ..Continuations: find_errors

function find_errors(solver::Codim2Solver)
    if augsys(solver.prob) isa BackReferencingAS
        # TODO: find_errors for BackReferencingAS-problems
        @warn("`find_errors` for problem with `BackReferencingAS`",
             " is not implemented yet.")
        return []
    end
    super = as(solver, ContinuationSolver)
    errors = find_errors(super)
    # TODO: do some Codim2Solver-specific checks
    return errors
end
