export RingPolymerPhasespace

"""
    RingPolymerPhasespace{T} <: RingPolymerDynamicalVariables{T}

Container for the dynamical positions and momenta of the atoms in the system.
"""
struct RingPolymerPhasespace{T} <: RingPolymerDynamicalVariables{T}
    x::Array{T, 3}
    RingPolymerPhasespace{T}(z::Array{N,3}) where {T<:AbstractFloat} where N = new(austrip.(z))
end

RingPolymerPhasespace(RP::Array{T, 3}) where {T<:Real} = RingPolymerPhasespace{T}(RP)

function RingPolymerPhasespace(R::Matrix{T}, P::Matrix{T}) where {T<:Real}
    RingPolymerPhasespace{T}(permutedims(cat(R, P, dims=3), [1, 3, 2]))
end

function RingPolymerPhasespace(
    R::Vector{T},
    P::Vector{T},
    beads::Integer,
    unit::Tuple{<:Unitful.Units,<:Unitful.Units}=(u"bohr", u"ħ_au/bohr")) where {T<:Real}
    RingPolymerPhasespace{T}(repeat(hcat(R.*unit[1], P.*unit[2]), 1, 1, beads))
end

get_positions(z::RingPolymerPhasespace) = @view z.x[:,1,:]
get_momenta(z::RingPolymerPhasespace) = @view z.x[:,2,:]