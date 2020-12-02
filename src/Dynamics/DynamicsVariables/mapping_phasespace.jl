using RecursiveArrayTools

export MappingPhasespace
export get_mapping_positions
export get_mapping_momenta

"""
    MappingPhasespace{T} <: DynamicalVariables{T}
    
Phasespace type that includes mapping variables.
"""
struct MappingPhasespace{T} <: DynamicalVariables{T}
    x::ArrayPartition{T}
end

"""
    MappingPhasespace(R::Vector{T}, P::Vector{T}, n_states::Integer, state::Integer, n_beads::Integer=1) where {T<:AbstractFloat}
    
Constructor for nonequilibrium mapping phasespace initialised on diabatic state `state`.
"""
function MappingPhasespace(R::Vector{T}, P::Vector{T}, n_states::Integer, state::Integer, n_beads::Integer=1) where {T<:AbstractFloat}
    z = Phasespace(R, P) 
    
    QP = zeros(T, n_states * n_beads, 2)
    θ = rand(T, n_states * n_beads) .* 2π
    QP[:,1] .= 1 ./ sqrt.(tan.(-θ) .^ 2 .+ 1)
    QP[state*(n_beads-1)+1:state*n_beads,1] .*= sqrt(3)
    QP[:,2] = QP[:,1] .* tan.(-θ)

    MappingPhasespace{T}(ArrayPartition(z.x, QP))
end

get_positions(z::MappingPhasespace) = @view z.x.x[1][:,1]
get_momenta(z::MappingPhasespace) = @view z.x.x[1][:,2]

get_mapping_positions(z::MappingPhasespace) = z.x.x[2][:,1]
# get_mapping_positions(z::MappingPhasespace, i::Integer) = z.x.x[2][:,1]

get_mapping_momenta(z::MappingPhasespace) = z.x.x[2][:,2]
# get_mapping_momenta(z::MappingPhasespace, i::Integer) = z.x.x[2][:,2]