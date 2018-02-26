struct ContinuationSweep{uType <: AbstractVector}
    u::Vector{uType}
    simple_bifurcation::Vector{Int}
end

ContinuationSweep(uType::Type{<: AbstractArray}) =
    ContinuationSweep(Vector{uType}(),
                      Vector{Int}())

struct ContinuationSolution{uType <: AbstractVector}
    sweeps::Vector{ContinuationSweep{uType}}
end

ContinuationSolution(uType::Type{<: AbstractArray}) =
    ContinuationSolution([ContinuationSweep(uType)])

function push_point!(sweep::ContinuationSweep, u,
                     simple_bifurcation::Bool = false)
    if simple_bifurcation
        push!(sweep.simple_bifurcation, length(sweep.u))
    end
    push!(sweep.u, u)
end

function push_point!(sol::ContinuationSolution, args...)
    push_point!(sol.sweeps[end], args...)
end

function new_sweep!(sol::ContinuationSolution{uType}) where uType
    push!(sol.sweeps, ContinuationSweep(uType))
end
