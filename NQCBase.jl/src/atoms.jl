
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
struct Atoms{Vt<:AbstractVector{Symbol}, Vn<:AbstractVector{Int}, Vm<:AbstractVector}
    types::Vt
    numbers::Vn
    masses::Vm
end

# Default to CPU-based Atoms
function Atoms(atom_types::AbstractVector{Symbol})
    Atoms(atom_types, CPUDevice)
end

# CPU-based Atoms constructor. Use Float64, Int64 by default.
function Atoms(atom_types::AbstractVector{Symbol}, CPUDevice)
    types = Vector{Symbol}(atom_types)
    numbers = Vector{Int}([element.number for element in elements[atom_types]])
    masses = Vector{Float64}([austrip(element.atomic_mass) for element in elements[atom_types]])
    return Atoms(types, numbers, masses)
end

# Default use of CPUs
Atoms(masses::AbstractVector) = Atoms(masses, CPUDevice)

function Atoms(masses::AbstractVector, CPUDevice)
    natoms = length(masses)
    types = Vector{Symbol}(fill(:X, natoms))
    numbers = Vector{Int}(zeros(natoms))
    masses = Vector{Float64}(austrip.(masses))
    Atoms(types, numbers, masses)
end

Atoms(atom_type) = Atoms([atom_type])

Base.length(atoms::Atoms) = length(atoms.types)
Base.range(atoms::Atoms) = range(1; length=length(atoms))

Base.IndexStyle(::Atoms) = IndexLinear()
Base.getindex(A::Atoms, i::Int) = Atoms(A.types[i])
Base.getindex(A::Atoms, i::AbstractRange) = Atoms(A.types[i])

function Base.:(==)(a::Atoms, b::Atoms)
    (a.types == b.types) && (a.numbers == b.numbers) && (a.masses ≈ b.masses)
end
