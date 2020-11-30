export DynamicalVariables
export Phasespace
export get_positions
export get_momenta

"""
Abstract type for different kinds of systems.
"""
DynamicalVariables{T} = DEDataVector{T} where T<:AbstractFloat

"""
    Phasespace{T} <: DynamicalVariables{T}

Container for the dynamical positions and momenta of the atoms in the system.
"""
struct Phasespace{T} <: DynamicalVariables{T}
    x::ArrayPartition{T, Tuple{Vector{T}, Vector{T}}}
    Phasespace{T}(z::ArrayPartition) where {T<:AbstractFloat} = new(austrip.(z))
end

function Phasespace(x::ArrayPartition{T, Tuple{Vector{T}, Vector{T}}}) where T<:AbstractFloat
    Phasespace{T}(x)
end

function Phasespace(
    R::Vector{T},
    P::Vector{T},
    unit::Tuple{<:Unitful.Units,<:Unitful.Units}=(u"bohr", u"ħ_au/bohr")) where {T<:Real}
    Phasespace{T}(ArrayPartition(R.*unit[1], P.*unit[2]))
end

get_positions(z::Phasespace) = z.x.x[1]
get_momenta(z::Phasespace) = z.x.x[2]
function Base.zero(z::Phasespace)
    blank = deepcopy(z)
    blank.x .= zero(z.x)
    blank
end