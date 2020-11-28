export Cell

"""
    Cell{T<:Unitful.Length{<:AbstractFloat}}

Simulation cell
"""
struct Cell{T<:AbstractFloat}
    vectors::Matrix{<:Unitful.Length{T}}
end

function Cell(vectors::Matrix{T}, unit::Unitful.Units=u"bohr") where {T<:AbstractFloat}
    Cell(vectors .* unit) 
end

function Cell(vectors::Matrix{T}, unit::Unitful.Units=u"bohr") where {T<:Integer}
    Cell(vectors .* unit .* 1.0) 
end
