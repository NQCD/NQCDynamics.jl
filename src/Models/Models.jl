module Models

using NeighbourLists
using StaticArrays
using UnitfulAtomic
using LinearAlgebra

using ..NonadiabaticMolecularDynamics

export Model
export AdiabaticModel
export DiabaticModel

abstract type Model end
abstract type AdiabaticModel <: Model end
abstract type DiabaticModel <: Model end

"""
    get_pairs(R::AbstractMatrix, cutoff::AbstractFloat, cell::AbstractCell)
    
Use NeighbourLists to calculate the neighbour list.
"""
function get_pairs(R::AbstractMatrix, cutoff::AbstractFloat, cell::AbstractCell)
    Q = copy(reinterpret(SVector{size(R)[1], eltype(R)}, vec(R))) # Convert array to vector of SVectors
    PairList(Q, cutoff, austrip.(cell.vectors'), cell.periodicity) # Construct the list of neighbours
end

include("analytic_models/free.jl")
include("analytic_models/harmonic.jl")
include("analytic_models/double_well.jl")
include("analytic_models/tully_models.jl")
include("EAM/PdH.jl")
include("EANN/EANN_H2Cu.jl")
include("EANN/EANN_H2Ag.jl")

include("plot.jl")
end # module
