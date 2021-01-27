export AbstractCell
export PeriodicCell
export InfiniteCell
export set_periodicity!
export apply_cell_boundaries!

abstract type AbstractCell end

struct InfiniteCell <: AbstractCell end

"""
    PeriodicCell{T<:AbstractFloat} <: AbstractCell

Optionally periodic cell
"""
struct PeriodicCell{T<:AbstractFloat} <: AbstractCell
    vectors::Matrix{T}
    inverse::Matrix{T}
    periodicity::Vector{Bool}
    function PeriodicCell{T}(vectors::Matrix, periodicity::Vector{Bool}) where {T}
        new(vectors, inv(vectors), periodicity)
    end
end

function PeriodicCell(vectors::AbstractMatrix{T}) where {T<:AbstractFloat}
    PeriodicCell{T}(vectors, [true, true, true]) 
end

function PeriodicCell(vectors::AbstractMatrix{<:Integer})
    PeriodicCell{Float64}(vectors, [true, true, true]) 
end

function set_periodicity!(cell::PeriodicCell, periodicity::Vector{Bool})
    cell.periodicity .= periodicity
end

function apply_cell_boundaries!(cell::PeriodicCell, R::AbstractMatrix)
    @views for i in axes(R, 2) # atoms
        R[:,i] .= cell.inverse * R[:,i]
        R[:,i] .= mod1.(R[:,i], 1)
        R[:,i] .= cell.vectors * R[:,i]
    end
end
apply_cell_boundaries!(::InfiniteCell, ::AbstractMatrix) = nothing

"""
    apply_cell_boundaries!(cell::PeriodicCell, R::AbstractArray{T,3}, beads::RingPolymerParameters) where {T}
    
Apply cell boundaries to the ring polymer system.

This converts to normal mode coordinates and restricts only the centroid.
This means that replicas can exit the cell but the centroid cannot.
The reasoning for this is to avoid complications in the computation of the spring potential.
Proper treatment would require accounting for the periodic cell in this function which I have
not yet done.

The modification by a factor of ``\\sqrt{N}`` is to convert to real space centroid. 
"""
function apply_cell_boundaries!(cell::PeriodicCell, R::AbstractArray{T,3}, beads::RingPolymerParameters) where {T}
    transform_to_normal_modes!(beads, R)
    R[:,beads.quantum_atoms,1] ./= sqrt(length(beads))
    @views apply_cell_boundaries!(cell, R[:,:,1])
    R[:,beads.quantum_atoms,1] .*= sqrt(length(beads))
    transform_from_normal_modes!(beads, R)
end
apply_cell_boundaries!(::InfiniteCell, ::AbstractArray{T,3}, ::RingPolymerParameters) where {T} = nothing