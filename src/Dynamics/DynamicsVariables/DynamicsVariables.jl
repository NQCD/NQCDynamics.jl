module DynamicsVariables

using RecursiveArrayTools
using Unitful
using UnitfulAtomic
using DiffEqBase

export DynamicalVariables
export RingPolymerDynamicalVariables
export get_positions
export get_momenta

"""
Abstract type for different kinds of systems.
"""
DynamicalVariables{T} = DEDataVector{T} where {T<:AbstractFloat}

abstract type RingPolymerDynamicalVariables{T} <: DynamicalVariables{T} end

include("phasespace.jl")
include("mapping_phasespace.jl")
include("surface_hopping_phasespace.jl")
include("ring_polymer_phasespace.jl")

end