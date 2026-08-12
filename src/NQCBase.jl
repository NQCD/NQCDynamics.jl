module NQCBase

using PeriodicTable
using Unitful, UnitfulAtomic
using Requires

# Basic types for atomic structures
include("unit_conversions.jl")
include("atoms.jl")
include("cells.jl")

# # Marker types for adiabatic and diabatic states.
# Used in NQCModels to mark which type of Hamiltonian is provided.
# Used in NQCDistributions to mark whether an initial electronic state is adiabatic or diabatic.
abstract type StateType end

struct Adiabatic <: StateType end
struct Diabatic <: StateType end

export Adiabatic, Diabatic

# Convenience Structure type

"""
    Structure{T}

Structures are storage types to keep atoms, positions, cells and further information in one place.
Many functions in the NQCD packages use only parts of the structure information for efficiency, e.g. integration routines, but it can be convenient to have it all in one type.

If you are developing in NQCD, please include multiple dispatch versions of your functions using Structure types where this would be convenient to the user.
"""
struct Structure
    atoms::NQCBase.Atoms # Atoms object
    positions::AbstractMatrix # Positions matrix in atomic units, each column containing one atom's positions
    cell::AbstractCell # Unit cell object
    info::Dict{String, Any} # Other structure information, e.g. from an ExtXYZ header
end
function Structure(atoms::NQCBase.Atoms, positions::AbstractMatrix, cell::AbstractCell)
    @assert length(atoms.types) == size(positions, 2) "Size of Positions needs to match number of Atoms. "
    if isa(cell, PeriodicCell)
        @assert size(cell.vectors, 1) == size(positions, 1) "Unit cell vectors must have the same dimensionality as Positions. "
    end
    return Structure(atoms, positions, cell, Dict{String, Any}())
end

# I/O Interfaces
include("io/extxyz.jl")
include("atoms_base.jl")


export Cell
export System
export Trajectory
export Position, Velocity





end
