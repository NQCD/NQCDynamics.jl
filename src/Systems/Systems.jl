module Systems

using RecursiveArrayTools
using PeriodicTable
using Unitful
using UnitfulAtomic
using DiffEqBase

include("cell.jl")
include("parameters.jl")
include("phasespace.jl")

export System

"""
    System{T<:AbstractFloat}

Simple container for both the static parameters and dynamical variables
"""
struct System{T<:AbstractFloat}
    parameters::SystemParameters{T}
    dynamical_variables::DynamicalVariables{T}
end

end # module


