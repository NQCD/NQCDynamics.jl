using LinearAlgebra: norm, mul!
using Distances: evaluate, PeriodicEuclidean

export AbstractCell
export PeriodicCell
export Supercell
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
    function PeriodicCell{T}(vectors::AbstractMatrix, periodicity::Vector{Bool}) where {T}
        new{T}(vectors, inv(vectors), periodicity,
            zeros(size(vectors)[1]), zeros(size(vectors)[1]), zeros(Bool, size(vectors)[1]))
    end
end

Base.eltype(::PeriodicCell{T}) where {T} = T

function PeriodicCell(vectors::AbstractMatrix)
    vectors = austrip.(vectors)
    PeriodicCell{eltype(vectors)}(vectors, [true, true, true]) 
end

function PeriodicCell(vectors::AbstractMatrix{<:Integer})
    PeriodicCell{Float64}(vectors, [true, true, true]) 
end

function set_periodicity!(cell::PeriodicCell, periodicity::AbstractVector{Bool})
    cell.periodicity .= periodicity
end

function set_vectors!(cell::PeriodicCell, vectors::AbstractMatrix)
    cell.vectors .= vectors
    cell.inverse .= inv(cell.vectors)
end

function apply_cell_boundaries!(cell::PeriodicCell, R::AbstractMatrix)
    @views for i in axes(R, 2) # atoms
        apply_cell_boundaries!(cell, R[:,i])
    end
end
apply_cell_boundaries!(::InfiniteCell, ::AbstractArray) = nothing

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

# Handler for periodic replicas of structures
struct PeriodicReplica{T}
    cell::PeriodicCell{T}
    translation::AbstractVector{Int}
end

function (operation::PeriodicReplica)(positions::AbstractMatrix)
    pos_translated = copy(positions)
    @inbounds for idx in axes(pos_translated, 2)
        pos_translated[:, idx] .+= operation.cell.vectors * operation.translation
    end
    return pos_translated
end

"""
    Supercell{T}(cell::PeriodicCell, replicas::AbstractVector)

This type carries information about how to replicate periodic copies of a structure based on the unit cell provided. 
"""
struct Supercell{T}
    cell::PeriodicCell{T}
    replicas::AbstractVector{PeriodicReplica{T}}
end

"""
    Supercell(cell::PeriodicCell, x_range, y_range, z_range)

Creates a `Supercell` which can be applied to a matrix of positions to generate all periodic copies, e.g. for plotting purposes in the following way:
```
sc = Supercell(cell, -1:1, [1,2], 1)
sc(positions) # Gives a total matrix containing the initial positions and all possible translations in the respective directions
```

Three iterators must be supplied, for which all possible translations will be generated. These can be integers (e.g. leave `z_range = 0` to generate no replicas in z-direction), 
or any other iterable that produces an Integer. 

All translations have to be (signed) integers, so [1,0,0] is one unit cell translation along the positive x direction. 

## Arguments

- `cell`: Specification for a unit cell
- `x_range`: Iterator yielding `Int`s to specify all x translations. 
- `y_range`: Iterator yielding `Int`s to specify all y translations. 
- `z_range`: Iterator yielding `Int`s to specify all z translations. 
"""
function Supercell(cell::PeriodicCell, x_range, y_range, z_range)
    # Build all possible combinations for replication
    replicas = [PeriodicReplica(cell, vcat(i...)) for i in Iterators.product(x_range, y_range, z_range)] |> vec
    return Supercell(cell, replicas)
end

function (supercell::Supercell)(positions::AbstractMatrix)
    return hcat(positions, [Rep(positions) for Rep in supercell.replicas]...)
end