using PeriodicTable
using UnitfulAtomic
using StaticArrays

export Atoms

struct Atoms{S,T<:AbstractFloat}
    types::SVector{S,Symbol}
    numbers::SVector{S,UInt8}
    masses::SVector{S,T}
end

function Atoms{T}(atom_types::AbstractVector{Symbol}) where {T}
    S = length(atom_types)
    types = SVector{S,Symbol}(atom_types)
    numbers = SVector{S,UInt8}([element.number for element in elements[atom_types]])
    masses = SVector{S,T}([austrip(element.atomic_mass) for element in elements[atom_types]])
    Atoms{S,T}(types, numbers, masses)
end

function Atoms{T}(masses::AbstractVector) where {T}
    S = length(masses)
    types = SVector{S,Symbol}(fill(:X, S))
    numbers = SVector{S,UInt8}(zeros(S))
    masses = SVector{S,T}(austrip.(masses))
    Atoms{S,T}(types, numbers, masses)
end

Atoms{T}(atom_type) where {T} = Atoms{T}([atom_type])
Atoms(atom_types) = Atoms{Float64}(atom_types)

Base.length(::Atoms{S}) where {S} = S
Base.range(::Atoms{S}) where {S} = range(1; length=S)

Base.IndexStyle(::Atoms) = IndexLinear()
Base.getindex(A::Atoms{S,T}, i::Int) where {S,T} = Atoms{T}(A.types[i])
Base.getindex(A::Atoms{S,T}, i::AbstractRange) where {S,T} = Atoms{T}(A.types[i])
