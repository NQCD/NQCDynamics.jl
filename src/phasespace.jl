
using RecursiveArrayTools
using DiffEqBase

export DynamicalVariables
export Phasespace
export RingPolymerPhasespace
export get_positions
export get_momenta
export get_flat_positions
export get_flat_momenta

"""
Abstract type for different kinds of systems.
"""
abstract type DynamicalVariables{T<:AbstractFloat} <: DEDataVector{T} end

mutable struct Phasespace{T,N} <: DynamicalVariables{T}
    x::ArrayPartition{T, Tuple{Array{T,N}, Array{T,N}}}
end
Phasespace(R::Matrix, P::Matrix) = Phasespace(ArrayPartition(R, P))
Phasespace(R::Matrix) = Phasespace(R, R)

const RingPolymerPhasespace{T} = Phasespace{T,3}
RingPolymerPhasespace(R::Array{T, 3}, P::Array{T, 3}) where {T} = Phasespace{T,3}(ArrayPartition(R, P))
function RingPolymerPhasespace(R::Matrix, P::Matrix, n_beads::Integer)
    R = cat([R for i=1:n_beads]..., dims=3)
    P = cat([P for i=1:n_beads]..., dims=3)
    RingPolymerPhasespace(R, P)
end

get_positions(z::DynamicalVariables) = z.x.x[1]
get_momenta(z::DynamicalVariables) = z.x.x[2]

get_flat_positions(z::Phasespace) = @view z[1:length(z)÷2]
get_flat_momenta(z::Phasespace) = @view z[length(z)÷2+1:end]

"""Get the positions for bead `i`"""
get_positions(z::RingPolymerPhasespace, i::Integer) = @view get_positions(z)[:,:,i]
get_momenta(z::RingPolymerPhasespace, i::Integer) = @view get_momenta(z)[:,:,i]