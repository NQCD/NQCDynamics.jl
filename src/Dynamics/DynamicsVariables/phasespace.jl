export Phasespace
export get_positions
export get_momenta

export RingPolymerPhasespace
export get_bead_positions
export get_bead_momenta

struct Phasespace{T} <: DynamicalVariables{T}
    x::ArrayPartition{T, Tuple{Matrix{T}, Matrix{T}}}
end
Phasespace(R::Matrix, P::Matrix) = Phasespace(ArrayPartition(R, P))

const RingPolymerArray{T, N} = ArrayPartition{T, NTuple{N, Matrix{T}}} where {T,N}

struct RingPolymerPhasespace{T,N} <: DynamicalVariables{T}
    x::ArrayPartition{T, Tuple{RingPolymerArray{T,N}, RingPolymerArray{T,N}}}
end
RingPolymerPhasespace(R::ArrayPartition, P::ArrayPartition) = RingPolymerPhasespace(ArrayPartition(R, P))
function RingPolymerPhasespace(R::Matrix, P::Matrix, n_beads::Integer)
    R = ArrayPartition([R for i=1:n_beads]...)
    P = ArrayPartition([P for i=1:n_beads]...)
    RingPolymerPhasespace(R, P)
end

get_positions(z::Union{RingPolymerPhasespace, Phasespace}) = z.x.x[1]
get_momenta(z::Union{RingPolymerPhasespace, Phasespace}) = z.x.x[2]

"""Get the positions for bead `i`"""
get_positions(z::RingPolymerPhasespace, i::Integer) = get_positions(z).x[i]
get_momenta(z::RingPolymerPhasespace, i::Integer) = get_momenta(z).x[i]

"""Get the postions of all the beads for DoF `i`"""
get_bead_positions(z::RingPolymerPhasespace, i::Integer, n_DoF::Integer) = @view get_positions(z)[i:n_DoF:end]
get_bead_momenta(z::RingPolymerPhasespace, i::Integer, n_DoF::Integer) = @view get_momenta(z)[i:n_DoF:end]