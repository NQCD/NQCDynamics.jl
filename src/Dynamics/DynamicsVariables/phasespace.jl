export Phasespace

"""
    Phasespace{T} <: DynamicalVariables{T}

Container for the dynamical positions and momenta of the atoms in the system.
"""
struct Phasespace{T} <: DynamicalVariables{T}
    x::Matrix{T}
    Phasespace{T}(z::Matrix) where {T<:AbstractFloat} = new(austrip.(z))
end

Phasespace(RP::Matrix{T}) where {T<:Real} = Phasespace{T}(RP)

function Phasespace(
    R::Vector{T},
    P::Vector{T},
    unit::Tuple{<:Unitful.Units,<:Unitful.Units}=(u"bohr", u"ħ_au/bohr")) where {T<:Real}
    Phasespace{T}(hcat(R.*unit[1], P.*unit[2]))
end

get_positions(z::Phasespace) = @view z.x[:,1]
get_momenta(z::Phasespace) = @view z.x[:,2]