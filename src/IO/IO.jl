module IO

using PyCall
using Unitful
using ..Atoms
using ..Systems
using ..Dynamics
using DiffEqBase

export read_system
export write_trajectory

function read_system(file::String)
    io = pyimport("ase.io")
    atoms = io.read(file)

    cell = Systems.PeriodicCell(atoms.cell.data, u"Å")
    atom_types = Symbol.(atoms.get_chemical_symbols())
    p = Atoms.AtomicParameters(cell, atom_types)

    R = vcat(atoms.get_positions()'...)
    P = zero(R)
    z = Dynamics.Phasespace(R, P)
    
    (p, z)
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