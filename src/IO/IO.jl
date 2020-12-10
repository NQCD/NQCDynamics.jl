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

function read_system(file::String)
    io = pyimport("ase.io")
    atoms = io.read(file)

    extract_parameters_and_positions(atoms)
end

function extract_parameters_and_positions(atoms::PyObject)
    cell = Systems.PeriodicCell(atoms.cell.data, u"Å")
    atom_types = Symbol.(atoms.get_chemical_symbols())
    p = Atoms.AtomicParameters(cell, atom_types)

    R = austrip.(atoms.get_positions()'u"Å")
    (p, R)
end

function write_trajectory(file_name::String, solution::DESolution, p::AtomicParameters)
    ase = pyimport("ase")
    trajectory = []
    for coord in get_positions.(solution.u)
        R = reshape(coord, 3, p.n_atoms)'
        atoms = ase.Atoms(positions=R)
        atoms.set_chemical_symbols(string.(p.atom_types))
        push!(trajectory, atoms)
    end

    ase.io.write(file_name, trajectory)
end

end # module