using LinearAlgebra: norm, mul!
using Distances: evaluate, PeriodicEuclidean

export AbstractCell
export PeriodicCell
export InfiniteCell
export set_periodicity!
export set_vectors!
export apply_cell_boundaries!
export evaluate_periodic_distance
export check_atoms_in_cell

const periodic_distance = PeriodicEuclidean([1, 1, 1])

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
    tmp_vector1::Vector{T}
    tmp_vector2::Vector{T}
    tmp_bools::Vector{Bool}
    function PeriodicCell{T}(vectors::Matrix, periodicity::Vector{Bool}) where {T}
        new(vectors, inv(vectors), periodicity,
            zeros(size(vectors)[1]), zeros(size(vectors)[1]), zeros(Bool, size(vectors)[1]))
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

function set_vectors!(cell::PeriodicCell, vectors::Matrix)
    cell.vectors .= vectors
    cell.inverse .= inv(cell)
end

function apply_cell_boundaries!(cell::PeriodicCell, R::AbstractMatrix)
    @views for i in axes(R, 2) # atoms
        mul!(cell.tmp_vector1, cell.inverse, R[:,i])
        for j in axes(R, 1) # DoFs
            if cell.periodicity[j]
                cell.tmp_vector1[j] = mod(cell.tmp_vector1[j], 1)
            end
        end
        mul!(R[:,i], cell.vectors, cell.tmp_vector1)
    end
end
apply_cell_boundaries!(::InfiniteCell, ::AbstractMatrix) = nothing

function apply_cell_boundaries!(cell::PeriodicCell, R::AbstractVector)
    mul!(cell.tmp_vector1, cell.inverse, R)
    for j in axes(R, 1) # DoFs
        if cell.periodicity[j]
            cell.tmp_vector1[j] = mod(cell.tmp_vector1[j], 1)
        end
    end
    mul!(R, cell.vectors, cell.tmp_vector1)
end

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

"""
    check_atoms_in_cell(cell::PeriodicCell, R::AbstractMatrix)::Bool

True if all atoms are inside the cell, false otherwise.
"""
function check_atoms_in_cell(cell::PeriodicCell, R::AbstractMatrix)::Bool
    @views for i in axes(R, 2) # atoms
        mul!(cell.tmp_vector1, cell.inverse, R[:,i])
        @. cell.tmp_bools = (cell.tmp_vector1 > 1) | (cell.tmp_vector1 < 0)
        any(cell.tmp_bools) && return false
    end
    true
end

function evaluate_periodic_distance(cell::PeriodicCell, r1::AbstractVector, r2::AbstractVector)
    mul!(cell.tmp_vector1, cell.inverse, r1)
    mul!(cell.tmp_vector2, cell.inverse, r2)
    evaluate(periodic_distance, cell.tmp_vector1, cell.tmp_vector2)
end