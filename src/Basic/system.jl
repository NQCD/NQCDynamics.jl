export SystemParameters

"""
    SystemParameters{T<:AbstractFloat}

Parameters for the system.

Contains all the static information about the system:
The simulation cell, the atom types and the number of atoms/ring polymer beads.
"""
struct SystemParameters{T<:AbstractFloat}
    cell::Cell{T}
    atom_types::Vector{Symbol}
    masses::Vector{T}
    n_atoms::UInt
    n_beads::UInt
end

"""
    SystemParameters(
        cell::Cell{<:Unitful.Length{T}},
        atom_types::Vector{Symbol}, 
        n_beads::Integer=1) where {T<:AbstractFloat}

Constructor for the parameters.

The user must provide the simulation cell and the list of elements in the system.
Optionally they can provide a number of beads for ring polymer simulations.
"""
function SystemParameters(
    cell::Cell{T},
    atom_types::Vector{Symbol}, 
    n_beads::Integer=1) where {T<:AbstractFloat}

    masses = austrip.([elements[atom].atomic_mass for atom in atom_types])
    SystemParameters{T}(cell, atom_types, masses, length(atom_types), n_beads)
end