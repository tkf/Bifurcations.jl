struct ContinuationSolution{U <: AbstractVector,
                            SB <: AbstractVector{<: Integer}}
    u::U
    simple_bifurcation::SB
end

ContinuationSolution(uType::Type{<: AbstractArray}) =
    ContinuationSolution(Vector{uType}(),
                         Vector{Int}())

function push_point!(sol::ContinuationSolution, u, simple_bifurcation::Bool)
    if simple_bifurcation
        push!(sol.simple_bifurcation, length(sol.u))
    end
    push!(sol.u, u)
end
