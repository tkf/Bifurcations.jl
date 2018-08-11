stabilities(sweep::Codim1Sweep) = isstable.((sweep.timekind,), sweep.eigvals)

regions(xs::AbstractVector{<: Bool}) = regions(identity, xs)

function regions(f, xs::AbstractVector)
    ranges = UnitRange{Int}[]
    if isempty(xs)
        return ranges
    end

    s = 1
    b = f(xs[1])
    @inbounds for i in 2:length(xs)
        c = f(xs[i])
        if c != b
            push!(ranges, s:i-1)
            s = i
            b = c
        end
    end
    if s !== length(xs)
        push!(ranges, s:length(xs))
    end
    return ranges
end

# TODO: Merge curves_by_stability to curves(...; by=:stability)
function curves_by_stability(sweep::Codim1Sweep, vars)
    ss = stabilities(sweep)
    ranges = regions(ss)
    info = [ss[first(r)] for r in ranges]

    super = as(sweep, ContinuationSweep)
    TE = eltype(eltype(super))
    curves = Tuple(Vector{TE}[] for _ in vars)
    for r in ranges
        if first(r) == 1
            ran = r
        else
            ran = first(r) - 1:last(r)
        end
        for (i, v) in enumerate(vars)
            push!(curves[i], [u[v] for u in @view super.u[ran]])
        end
    end

    return info, curves
end

function curves_by_stability(sol::Codim1Solution)
    if isempty(sol.sweeps)
        return ([], [])
    end

    info, curves = curves_by_stability(sol.sweeps[1])
    for sweep in sol.sweeps[2:end]
        info2, curves2 = curves_by_stability(sweep)
        append!(info, info2)
        append!(curves, curves2)
    end
    return info, curves
end
