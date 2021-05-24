export EhrenfestVariables
export get_density_matrix

mutable struct EhrenfestVariables{T,D}  <: DynamicalVariables{T}
    x::ArrayPartition{Complex{T}, Tuple{D,D,Matrix{Complex{T}}}}
end

function EhrenfestVariables(x::ArrayPartition{T}) where {T<:AbstractFloat}
    EhrenfestVariables(ArrayPartition(x.x[1], x.x[2], Complex.(x.x[3])))
end

function EhrenfestVariables(v::AbstractArray, r::AbstractArray, n_states::Integer, state::Integer)
    σ = zeros(Complex{eltype(r)}, n_states, n_states)
    σ[state, state] = 1
    EhrenfestVariables(ArrayPartition(v, r, σ))
end

get_density_matrix(u::EhrenfestVariables) = u.x.x[3]
