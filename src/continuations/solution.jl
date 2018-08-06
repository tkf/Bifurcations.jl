struct SimpleBifurcationInterval{uType}
    i::Int
    u0::uType
    u1::uType
    h::Float64  # hType?
    direction::Int
    sweep::WeakRef
end


struct ContinuationSweep{uType <: AbstractVector}
    i::Int
    direction::Int
    u::Vector{uType}
    simple_bifurcation::Vector{SimpleBifurcationInterval{uType}}
    sol::WeakRef
end

ContinuationSweep(uType::Type{<: AbstractArray}, i, direction, sol) =
    ContinuationSweep(i,
                      direction,
                      Vector{uType}(),
                      Vector{SimpleBifurcationInterval{uType}}(),
                      sol)

Base.eltype(::Type{<: ContinuationSweep{uType}}) where uType = uType
Base.length(sweep::ContinuationSweep) = length(sweep.u)

struct ContinuationSolution{uType <: AbstractVector}
    sweeps::Vector{ContinuationSweep{uType}}
    prob
end

ContinuationSolution(uType::Type{<: AbstractArray}, prob) =
    ContinuationSolution(ContinuationSweep{uType}[], prob)

function push_point!(sweep::ContinuationSweep, cache)
    push_point!(sweep, cache.u)
    if cache.simple_bifurcation
        u0 = sweep.u[end-1]
        u1 = sweep.u[end]
        push!(sweep.simple_bifurcation, SimpleBifurcationInterval(
            length(sweep.u) - 1,
            u0,
            u1,
            cache.h,
            -cache.direction,  # re-flip to save direction at u0
            WeakRef(sweep),
        ))
    end
end

function push_point!(sweep::ContinuationSweep, u::AbstractArray)
    push!(sweep.u, copy(u))
end

function push_point!(sol::ContinuationSolution, args...)
    push_point!(sol.sweeps[end], args...)
end

function new_sweep!(sol::ContinuationSolution{uType},
                    direction) where uType
    push!(sol.sweeps, ContinuationSweep(
        uType,
        length(sol.sweeps) + 1,
        direction,
        WeakRef(sol)))
end

function sweeps_as_vectors(sol::ContinuationSolution{uType}, i) where uType
    if i == 0
        i = length(sol.sweeps[1].u[1])
    end
    x1 = [x[i] for x in sol.sweeps[1].u]
    x2 = [x[i] for x in sol.sweeps[2].u]
    curves = [vcat(reverse!(x2), x1)]
    for sweep in sol.sweeps[3:end]
        push!(curves, [x[i] for x in sweep.u])
    end
    return curves
end
