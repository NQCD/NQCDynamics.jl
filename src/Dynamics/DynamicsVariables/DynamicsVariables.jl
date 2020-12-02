module DynamicsVariables

using RecursiveArrayTools
using Unitful
using UnitfulAtomic
using DiffEqBase

export DynamicalVariables

"""
Abstract type for different kinds of systems.
"""
DynamicalVariables{T} = DEDataVector{T} where {T<:AbstractFloat}

include("phasespace.jl")
include("mapping_phasespace.jl")
include("surface_hopping_phasespace.jl")

end