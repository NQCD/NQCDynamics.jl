export AbstractCell
export PeriodicCell
export InfiniteCell
export set_periodicity!

abstract type AbstractCell{T<:AbstractFloat} end

"""
    struct Cell{T<:AbstractFloat} <: AbstractCell

Simulation cell
"""
struct PeriodicCell{T} <: AbstractCell{T}
    vectors::Matrix{<:Unitful.Length{T}}
    periodicity::Vector{Bool}
end

function PeriodicCell(vectors::Matrix{T}, unit::Unitful.Units=u"bohr") where {T<:AbstractFloat}
    PeriodicCell(vectors .* unit, [true, true, true]) 
end

function PeriodicCell(vectors::Matrix{T}, unit::Unitful.Units=u"bohr") where {T<:Integer}
    PeriodicCell(vectors .* unit .* 1.0, [true, true, true]) 
end

function set_periodicity!(cell::PeriodicCell, periodicity::Vector{Bool})
    cell.periodicity .= periodicity
end

struct InfiniteCell{T} <: AbstractCell{T} end
InfiniteCell() = InfiniteCell{Float64}()
