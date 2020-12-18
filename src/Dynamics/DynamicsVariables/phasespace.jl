export Phasespace
export get_positions
export get_momenta

export RingPolymerPhasespace

mutable struct Phasespace{T} <: DynamicalVariables{T}
    x::ArrayPartition{T, Tuple{Matrix{T}, Matrix{T}}}
end
Phasespace(R::Matrix, P::Matrix) = Phasespace(ArrayPartition(R, P))

mutable struct RingPolymerPhasespace{T} <: DynamicalVariables{T}
    x::ArrayPartition{T, Tuple{Array{T,3}, Array{T,3}}}
end
RingPolymerPhasespace(R::Array{T, 3}, P::Array{T, 3}) where {T} = RingPolymerPhasespace(ArrayPartition(R, P))
function RingPolymerPhasespace(R::Matrix, P::Matrix, n_beads::Integer)
    R = cat([R for i=1:n_beads]..., dims=3)
    P = cat([P for i=1:n_beads]..., dims=3)
    RingPolymerPhasespace(R, P)
end

get_positions(z::Union{RingPolymerPhasespace, Phasespace}) = z.x.x[1]
get_momenta(z::Union{RingPolymerPhasespace, Phasespace}) = z.x.x[2]

"""Get the positions for bead `i`"""
get_positions(z::RingPolymerPhasespace, i::Integer) = @view get_positions(z)[:,:,i]
get_momenta(z::RingPolymerPhasespace, i::Integer) = @view get_momenta(z)[:,:,i]