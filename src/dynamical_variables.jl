
using RecursiveArrayTools
using DiffEqBase

export DynamicalVariables
export get_positions
export get_velocities
export get_flat_positions
export get_flat_velocities
export ClassicalDynamicals

"""
Abstract type for different kinds of systems.
"""
abstract type DynamicalVariables{T<:AbstractFloat} <: DEDataVector{T} end

mutable struct ClassicalDynamicals{T,D} <: DynamicalVariables{T}
    x::ArrayPartition{T, Tuple{D, D}}
end
ClassicalDynamicals(v::AbstractArray, r::AbstractArray) = ClassicalDynamicals(ArrayPartition(v, r))
ClassicalDynamicals(v::AbstractArray) = ClassicalDynamicals(v, v)

get_velocities(z::DynamicalVariables) = z.x.x[1]
get_positions(z::DynamicalVariables) = z.x.x[2]

get_velocities(z::ArrayPartition) = z.x[1]
get_positions(z::ArrayPartition) = z.x[2]

get_flat_velocities(z::ClassicalDynamicals) = @view z[1:length(z)÷2]
get_flat_positions(z::ClassicalDynamicals) = @view z[length(z)÷2+1:end]
