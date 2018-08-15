import ..Codim1LimitCycle: limitcycles
using ..Codim1LimitCycle: AbstractLimitCycleData


function dim_ds_state(sweep::Codim2LCSweep)
    prob = as(sweep, ContinuationSweep).sol.value.prob
    return size(prob.xs0, 1)
end


struct FoldLimitCycleData{M, R, V} <: AbstractLimitCycleData
    state::M
    period::R
    param1_value::V
    param2_value::V
    # param1_axis::L
    # param2_axis::L
end

function FoldLimitCycleData(u::AbstractVector, n::Integer)
    state_elements, r = divrem(length(u) - 4, 2)
    @assert r == 0
    FoldLimitCycleData(
        reshape((@view u[1:state_elements]), n, :),
        u[state_elements + 1],  # period
        u[end - 1],             # param1_value
        u[end],                 # param2_value
    )
end
# See: [[./problem.jl::get_u0.*::]]

function limitcycles(sweep::Codim2LCSweep)
    n = dim_ds_state(sweep)
    u_list = as(sweep, ContinuationSweep).u
    return FoldLimitCycleData.(u_list, (n,))
end

limitcycles(sol::Codim2LCSolution) = vcat(limitcycles.(sol.sweeps)...)
limitcycles(solver::Codim2LCSolver) = limitcycles(solver.sol)


BifurcationsBase.measure(lc::FoldLimitCycleData, i::Integer) =
    @view lc.state[i, :]

# TODO: make it more type-based
function BifurcationsBase.measure(lc::FoldLimitCycleData, key::Symbol)
    if key == :p1
        return lc.param1_value
    elseif key == :p2
        return lc.param2_value
    elseif key == :period
        return lc.period
    end
    error("Unsupported key: $key")
end


# TODO: make it more type-based
function BifurcationsBase.measure(sweep::Codim2LCSweep, key::Symbol)
    u_list = as(sweep, ContinuationSweep).u
    if isempty(u_list)
        state_elements = 0  # won't be used
    else
        state_elements, r = divrem(length(u_list[1]) - 4, 2)
        @assert r == 0
    end
    if key == :p1
        return [u[end - 1] for u in u_list]
    elseif key == :p2
        return [u[end] for u in u_list]
    elseif key == :period
        return [u[state_elements + 1] for u in u_list]
    end
    error("Unsupported key: $key")
end
