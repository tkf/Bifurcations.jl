struct ContinuationSweep{uType <: AbstractVector}
    direction::Int
    u::Vector{uType}
    simple_bifurcation::Vector{Int}
end

ContinuationSweep(uType::Type{<: AbstractArray}, direction) =
    ContinuationSweep(direction,
                      Vector{uType}(),
                      Vector{Int}())

struct ContinuationSolution{uType <: AbstractVector}
    sweeps::Vector{ContinuationSweep{uType}}
end

ContinuationSolution(uType::Type{<: AbstractArray}) =
    ContinuationSolution(ContinuationSweep{uType}[])

function push_point!(sweep::ContinuationSweep, u,
                     simple_bifurcation::Bool = false)
    if simple_bifurcation
        push!(sweep.simple_bifurcation, length(sweep.u))
    end
    push!(sweep.u, copy(u))
end

function push_point!(sol::ContinuationSolution, args...)
    push_point!(sol.sweeps[end], args...)
end

function new_sweep!(sol::ContinuationSolution{uType},
                    direction) where uType
    push!(sol.sweeps, ContinuationSweep(uType, direction))
end

function sweeps_as_vectors(sol::ContinuationSolution{uType}, i) where uType
    if i == 0
        i = length(sol.sweeps[1].u[1])
    end
    x1 = [x[i] for x in sol.sweeps[1].u]
    x2 = [x[i] for x in sol.sweeps[2].u]
    return [vcat(reverse!(x2), x1)]
end
