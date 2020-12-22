module Calculators

using LinearAlgebra
using Unitful
using FiniteDifferences
using NeighbourLists
using StaticArrays
using UnitfulAtomic
using ..NonadiabaticMolecularDynamics: PeriodicCell

abstract type AbstractCalculator end

abstract type AdiabaticCalculator <: AbstractCalculator end

export evaluate_potential
export evaluate_derivative

function evaluate_potential(calc::AbstractCalculator, R::AbstractMatrix{T})::T where {T}
    calc.potential(R)
end

function evaluate_derivative(calc::AbstractCalculator, R::AbstractMatrix{T})::Matrix{T} where {T}
    calc.derivative(R)
end

"""
    get_pairs(R::AbstractMatrix, cutoff::AbstractFloat, cell::PeriodicCell)
    
Use NeighbourLists to calculate the neighbour list.
"""
function get_pairs(R::AbstractMatrix, cutoff::AbstractFloat, cell::PeriodicCell)
    Q = copy(reinterpret(SVector{size(R)[1], eltype(R)}, vec(R))) # Convert array to vector of SVectors
    PairList(Q, cutoff, cell.vectors', cell.periodicity) # Construct the list of neighbours
end

include("EAM/PdH.jl")
include("analytic/free.jl")
include("analytic/harmonic.jl")

end # module