export Phasespace
export get_positions
export get_momenta

"""
Abstract type for different kinds of systems.
    
We might also want to add alternatives to the phasespace that also store electronic
variables. An option to consider is the DEDataArray if we decide to go with 
DifferentialEquation.jl.
"""
abstract type DynamicalVariables{T<:AbstractFloat} end

"""
    Phasespace{T<:AbstractFloat} <: DynamicalVariables{T}

Container for the dynamical positions and momenta of the atoms in the system.
"""
struct Phasespace{T<:AbstractFloat} <: DynamicalVariables{T}
    z::ArrayPartition{T, Tuple{Vector{T}, Vector{T}}}
    Phasespace{T}(z::ArrayPartition) where {T<:AbstractFloat} = new(austrip.(z))
end

function Phasespace(
    R::Vector{T},
    P::Vector{T},
    unit::Tuple{<:Unitful.Units,<:Unitful.Units}=(u"bohr", u"ħ_au/bohr")) where {T<:Real}
    Phasespace{T}(ArrayPartition(R.*unit[1], P.*unit[2]))
end

get_positions(z::Phasespace) = z.z.x[1]
get_momenta(z::Phasespace) = z.z.x[2]