export AbstractCell
export PeriodicCell
export InfiniteCell
export set_periodicity!

abstract type AbstractCell{T<:AbstractFloat} end

struct InfiniteCell{T} <: AbstractCell{T} end

"""
    PeriodicCell{T} <: AbstractCell{T}

Optionally periodic cell
"""
struct PeriodicCell{T} <: AbstractCell{T}
    vectors::Matrix{T}
    periodicity::Vector{Bool}
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
