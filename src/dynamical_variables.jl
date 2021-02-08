
using RecursiveArrayTools
using DiffEqBase

export DynamicalVariables
export RingPolymerClassicalDynamicals
export get_positions
export get_velocities
export get_flat_positions
export get_flat_velocities
export ClassicalDynamicals

"""
Abstract type for different kinds of systems.
"""
abstract type DynamicalVariables{T<:AbstractFloat} <: DEDataVector{T} end

mutable struct ClassicalDynamicals{T} <: DynamicalVariables{T}
    x::ArrayPartition{T, Tuple{Matrix{T}, Matrix{T}}}
end
ClassicalDynamicals(v::Matrix, r::Matrix) = ClassicalDynamicals(ArrayPartition(v, r))
ClassicalDynamicals(v::Matrix) = ClassicalDynamicals(v, v)

abstract type RingPolymerDynamicalVariables{T} <: DynamicalVariables{T} end

mutable struct RingPolymerClassicalDynamicals{T} <: RingPolymerDynamicalVariables{T}
    x::ArrayPartition{T, Tuple{Array{T,3}, Array{T,3}}}
end
RingPolymerClassicalDynamicals(v::Array{T,3}, r::Array{T,3}) where {T} =
    RingPolymerClassicalDynamicals(ArrayPartition(v, r))

RingPolymerClassicalDynamicals(v::Array{T,3}) where {T} =
    RingPolymerClassicalDynamicals(v, v)

function RingPolymerClassicalDynamicals(v::Matrix, r::Matrix, n_beads::Integer)
    v = cat([v for i=1:n_beads]..., dims=3)
    r = cat([r for i=1:n_beads]..., dims=3)
    RingPolymerClassicalDynamicals(v, r)
end

get_velocities(z::DynamicalVariables) = z.x.x[1]
get_positions(z::DynamicalVariables) = z.x.x[2]

get_velocities(z::ArrayPartition) = z.x[1]
get_positions(z::ArrayPartition) = z.x[2]

get_flat_velocities(z::ClassicalDynamicals) = @view z[1:length(z)÷2]
get_flat_positions(z::ClassicalDynamicals) = @view z[length(z)÷2+1:end]

"""Get the positions for bead `i`"""
get_positions(z::RingPolymerClassicalDynamicals, i::Integer) = @view get_positions(z)[:,:,i]
get_velocities(z::RingPolymerClassicalDynamicals, i::Integer) = @view get_velocities(z)[:,:,i]