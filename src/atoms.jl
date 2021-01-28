using PeriodicTable
using UnitfulAtomic
using StaticArrays

export Atoms

struct Atoms{S,T<:AbstractFloat}
    types::SVector{S,Symbol}
    numbers::SVector{S,UInt8}
    masses::SVector{S,T}
    function Atoms{T}(atom_types::Vector{Symbol}) where {T}
        S = length(atom_types)
        types = SVector{S,Symbol}(atom_types)
        numbers = SVector{S,UInt8}([element.number for element in elements[atom_types]])
        masses = SVector{S,T}([austrip(element.atomic_mass) for element in elements[atom_types]])
        new{S,T}(types, numbers, masses)
    end
end

Atoms(atom_types::Vector{Symbol}) = Atoms{Float64}(atom_types)

Base.length(::Atoms{S}) where {S} = S
Base.range(::Atoms{S}) where {S} = range(1; length=S)