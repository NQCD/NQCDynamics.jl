module Models

using NeighbourLists
using StaticArrays
using UnitfulAtomic
using ..Atoms

export Model

abstract type Model end

"""
    get_pairs(R::Matrix, cutoff::AbstractFloat, cell::Atoms.AbstractCell)
    
Use NeighbourLists to calculate the neighbour list.
"""
function get_pairs(R::Matrix, cutoff::AbstractFloat, cell::Atoms.AbstractCell)
    Q = copy(reinterpret(SVector{size(R)[1], eltype(R)}, vec(R))) # Convert array to vector of SVectors
    PairList(Q, cutoff, austrip.(cell.vectors'), cell.periodicity) # Construct the list of neighbours
end

include("Analytic/Analytic.jl")
include("ML/ML.jl")
include("EAM/PdH.jl")
end # module
