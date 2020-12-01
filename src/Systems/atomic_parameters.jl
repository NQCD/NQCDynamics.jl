export AtomicParameters

"""
    AtomicParameters{T<:AbstractFloat}

Parameters for the system.

Contains all the static information about the system:
The simulation cell, the atom types and the number of atoms/ring polymer beads.
"""
struct AtomicParameters{T<:AbstractFloat}
    cell::AbstractCell{T}
    atom_types::Vector{Symbol}
    masses::Vector{T}
    n_atoms::UInt
end

"""
    AtomicParameters(cell::AbstractCell{T}, atom_types::Vector{Symbol}) where {T<:AbstractFloat}

Constructor for the parameters.

The user must provide the simulation cell and the list of elements in the system.
"""
function AtomicParameters(cell::AbstractCell{T}, atom_types::Vector{Symbol}) where {T<:AbstractFloat}

    masses = austrip.([elements[atom].atomic_mass for atom in atom_types])
    AtomicParameters{T}(cell, atom_types, masses, length(atom_types))
end