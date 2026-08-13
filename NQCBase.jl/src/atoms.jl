
export Atoms

"""
    Atoms{T<:AbstractFloat}

Basic atomic parameters: element symbols, numbers and masses

Masses are converted to atomic units.
Constructed using either element symbols or masses.

```jldoctest
julia> Atoms(:H)
Atoms{Float64}([:H], [1], [1837.4715941070515])

julia> Atoms([:H, :H, :H, :C])
Atoms{Float64}([:H, :H, :H, :C], [1, 1, 1, 6], [1837.4715941070515, 1837.4715941070515, 1837.4715941070515, 21894.713607956142])

julia> Atoms([100, 200])
Atoms{Float64}([:X, :X], [0, 0], [100.0, 200.0])
```
"""
struct Atoms{T<:AbstractFloat}
    types::Vector{Symbol}
    numbers::Vector{Int}
    masses::Vector{T}
end

function Atoms{T}(atom_types::AbstractVector{Symbol}) where {T}
    types = Vector{Symbol}(atom_types)
    numbers = Vector{Int}([element.number for element in elements[atom_types]])
    masses = Vector{T}([austrip(element.atomic_mass) for element in elements[atom_types]])
    Atoms{T}(types, numbers, masses)
end

function Atoms{T}(masses::AbstractVector) where {T}
    natoms = length(masses)
    types = Vector{Symbol}(fill(:X, natoms))
    numbers = Vector{Int}(zeros(natoms))
    masses = Vector{T}(austrip.(masses))
    Atoms{T}(types, numbers, masses)
end

Atoms{T}(atom_type) where {T} = Atoms{T}([atom_type])
Atoms(atom_types) = Atoms{Float64}(atom_types)

Base.length(atoms::Atoms) = length(atoms.types)
Base.range(atoms::Atoms) = range(1; length=length(atoms))

Base.IndexStyle(::Atoms) = IndexLinear()
Base.getindex(A::Atoms{T}, i::Int) where {T} = Atoms{T}(A.types[i])
Base.getindex(A::Atoms{T}, i::AbstractRange) where {T} = Atoms{T}(A.types[i])

function Base.:(==)(a::Atoms, b::Atoms)
    (a.types == b.types) && (a.numbers == b.numbers) && (a.masses ≈ b.masses)
end
