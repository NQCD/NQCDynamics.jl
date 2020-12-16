module IO

using PyCall
using Unitful
using UnitfulAtomic
using ..Atoms
using ..Systems
using ..Dynamics
using DiffEqBase

export read_system
export write_trajectory
export extract_parameters_and_positions
export create_ase_atoms

function read_system(file::String)
    io = pyimport("ase.io")
    atoms = io.read(file)

    extract_parameters_and_positions(atoms)
end

function extract_parameters_and_positions(atoms::PyObject)
    cell = Systems.PeriodicCell(atoms.cell.data', u"Å")
    atom_types = Symbol.(atoms.get_chemical_symbols())
    p = Atoms.AtomicParameters(cell, atom_types)

    R = austrip.(atoms.get_positions()'u"Å")
    (p, R)
end

function create_ase_atoms(atoms::AtomicParameters, R::Matrix{<:AbstractFloat}) 
    ase = pyimport("ase")
    if hasproperty(atoms.cell, :periodicity)
        return ase.Atoms(
            positions=ustrip.(u"Å", R'u"bohr"),
            cell=ustrip.(u"Å", atoms.cell.vectors),
            symbols=string.(atoms.atom_types),
            pbc=atoms.cell.periodicity)
    else
        return ase.Atoms(
            positions=ustrip.(u"Å", R'u"bohr"),
            symbols=string.(atoms.atom_types))
    end
end

function write_trajectory(file_name::String, solution::DESolution, p::AtomicParameters)
    ase = pyimport("ase")
    trajectory = []
    for coord in get_positions.(solution.u)
        R = reshape(coord, 3, p.n_atoms)
        atoms == create_ase_atoms(p, R)
        push!(trajectory, atoms)
    end

    ase.io.write(file_name, trajectory)
end

function write_trajectory(file_name::String, positions::Vector{Matrix{T}}, p::AtomicParameters) where {T}
    ase = pyimport("ase")
    trajectory = create_ase_atoms.(Ref(p), positions)

    ase.io.write(file_name, trajectory)
end

function write_trajectory(file_name::String, positions::Vector{Array{T, 3}}, p::AtomicParameters) where {T}
    ase = pyimport("ase")
    trajectory = []
    for config in positions
        push!(trajectory, sum([create_ase_atoms(p, config[:,:,i]) for i=1:size(positions[1])[3]]))
    end

    ase.io.write(file_name, trajectory)
end

end # module