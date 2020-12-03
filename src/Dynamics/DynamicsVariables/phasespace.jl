export Phasespace
export get_positions
export get_momenta
export get_bead_positions
export get_bead_momenta

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

function get_bead_positions(z::Phasespace, n_beads::Integer)
    collect(Iterators.partition(get_positions(z), n_beads))
end

get_momenta(z::Phasespace) = @view z.x[:,2]

function get_bead_momenta(z::Phasespace, n_beads::Integer)
    collect(Iterators.partition(get_momenta(z), n_beads))
end