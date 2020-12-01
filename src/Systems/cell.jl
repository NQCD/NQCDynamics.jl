export AbstractCell
export PeriodicCell
export InfiniteCell

abstract type AbstractCell{T<:AbstractFloat} end

"""
    struct Cell{T<:AbstractFloat} <: AbstractCell

Simulation cell
"""
struct PeriodicCell{T} <: AbstractCell{T}
    vectors::Matrix{<:Unitful.Length{T}}
end

function PeriodicCell(vectors::Matrix{T}, unit::Unitful.Units=u"bohr") where {T<:AbstractFloat}
    PeriodicCell(vectors .* unit) 
end

function PeriodicCell(vectors::Matrix{T}, unit::Unitful.Units=u"bohr") where {T<:Integer}
    PeriodicCell(vectors .* unit .* 1.0) 
end

struct InfiniteCell{T} <: AbstractCell{T} end
InfiniteCell() = InfiniteCell{Float64}()
