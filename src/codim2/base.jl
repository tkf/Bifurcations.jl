# ds_dim(x::AbstractVector) = ds_dim(length(x))
# ds_dim(var_dim) = var_dim ÷ 2 - 1
eq_dim(N) = 2N + 1
var_dim(N) = 2N + 2

function ds_state(x::AbstractArray)
    N = length(x) ÷ 2 - 1
    return @view x[1:N]
end

function ds_eigvec(x::AbstractArray)
    N = length(x) ÷ 2 - 1
    return @view x[N+1:2N]
end

function sarray_slice_epxr(C, x, indices)
    values = [:($x[$i]) for i in indices]
    return Expr(:call, C, values...)
end

@generated function ds_state(x::SVector{S, T}) where {S, T}
    N = S ÷ 2 - 1
    @assert var_dim(N) == S
    C = :(SVector{$N, T})
    return sarray_slice_epxr(C, :x, 1:N)
end

@generated function ds_eigvec(x::SVector{S, T}) where {S, T}
    N = S ÷ 2 - 1
    @assert var_dim(N) == S
    C = :(SVector{$N, T})
    return sarray_slice_epxr(C, :x, N+1:2N)
end

function output_vars(H)
    N = length(H) ÷ 2
    H1 = @view H[1:N]
    H2 = @view H[N+1:2N]
    H3 = @view H[end]
    return H1, H2, H3
end

cat_outputs(H1, H2, H3) = vcat(H1, H2, H3)

function cat_outputs(H1::SVector{S, T1},
                     H2::SVector{S, T2},
                     H3::T3) where {S, T1, T2, T3}
    values = [
        [:(H1[$i]) for i in 1:S]
        [:(H2[$i]) for i in 1:S]
        :(H3)
    ]
    T = promote_type(T1, T2, T3)
    return Expr(:call :(SVector{eq_dim(S), T}), values...)
end
