module Models

using NeighbourLists
using StaticArrays
using UnitfulAtomic
using LinearAlgebra
using Reexport

using ..NonadiabaticMolecularDynamics

export Model
export AdiabaticModel
export DiabaticModel
export FrictionModel

export potential!
export derivative!
export friction!

abstract type Model end
abstract type AdiabaticModel <: Model end
abstract type DiabaticModel <: Model end
abstract type FrictionModel <: Model end

"""
    get_pairs(R::AbstractMatrix, cutoff::AbstractFloat, cell::AbstractCell)
    
Use NeighbourLists to calculate the neighbour list.
"""
function get_pairs(R::AbstractMatrix, cutoff::AbstractFloat, cell::AbstractCell)
    Q = copy(reinterpret(SVector{size(R)[1], eltype(R)}, vec(R))) # Convert array to vector of SVectors
    PairList(Q, cutoff, austrip.(cell.vectors'), cell.periodicity) # Construct the list of neighbours
end

function potential! end
function derivative! end
function friction! end

include("analytic_models/free.jl")
include("analytic_models/harmonic.jl")

include("analytic_models/double_well.jl")
include("analytic_models/tully_models.jl")
include("analytic_models/scattering_anderson_holstein.jl")

include("analytic_models/friction_harmonic.jl")

include("EAM/PdH.jl")
include("EANN/EANN_H2Cu.jl")
include("EANN/EANN_H2Ag.jl")
include("ML/ML.jl")
@reexport using .ML

include("plot.jl")
end # module
