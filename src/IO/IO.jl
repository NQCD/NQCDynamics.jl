module IO

using PyCall
using Unitful
using UnitfulAtomic
using ..NonadiabaticMolecularDynamics
# using ..Atoms
# using ..Systems
# using ..Dynamics
# using DiffEqBase

export read_system
export write_trajectory
export extract_parameters_and_positions
export create_ase_atoms

function read_system(file::String)
    io = pyimport("ase.io")
    atoms = io.read(file)

    extract_parameters_and_positions(atoms)
end

function extract_parameters_and_positions(ase_atoms::PyObject)
    cell = PeriodicCell{Float64}(austrip.(ase_atoms.cell.data'u"Å"), ase_atoms.pbc)
    atom_types = Symbol.(ase_atoms.get_chemical_symbols())
    positions = austrip.(ase_atoms.get_positions()'u"Å")

    atoms = Atoms{Float64}(atom_types)
    (cell, atoms, positions)
end

function create_ase_atoms(cell::PeriodicCell, atoms::Atoms, R::Matrix) 
    ase = pyimport("ase")
    ase.Atoms(
        positions=ustrip.(u"Å", R'u"bohr"),
        cell=ustrip.(u"Å", cell.vectors'u"bohr"),
        symbols=string.(atoms.types),
        pbc=cell.periodicity)
end

create_ase_atoms(::InfiniteCell, atoms::Atoms, R::Matrix) = create_ase_atoms(atoms, R)

function create_ase_atoms(atoms::Atoms, R::Matrix) 
    ase = pyimport("ase")
    ase.Atoms(
        positions=ustrip.(u"Å", R'u"bohr"),
        symbols=string.(atoms.types))
end

# function write_trajectory(file_name::String, solution::DESolution, p::Atoms)
#     ase = pyimport("ase")
#     trajectory = []
#     for coord in get_positions.(solution.u)
#         R = reshape(coord, 3, p.n_atoms)
#         atoms == create_ase_atoms(p, R)
#         push!(trajectory, atoms)
#     end

#     ase.io.write(file_name, trajectory)
# end

function write_trajectory(file_name::String, cell::AbstractCell, atoms::Atoms, positions::Vector{<:Matrix})
    ase = pyimport("ase")
    trajectory = create_ase_atoms.(Ref(cell), Ref(atoms), positions)

    ase.io.write(file_name, trajectory)
end

function write_trajectory(file_name::String, cell::AbstractCell, atoms::Atoms, positions::Vector{Array{T, 3}}) where {T}
    trajectory = []
    for config in positions
        push!(trajectory, sum([create_ase_atoms(cell, atoms, config[:,:,i]) for i=1:size(positions[1])[3]]))
    end

    io = pyimport("ase.io")
    io.write(file_name, trajectory)
end

# function write_trajectory(file_name::String, positions::Vector{RingPolymerPhasespace}, p::Atoms) where {T}
#     ase = pyimport("ase")
#     trajectory = []
#     for config in positions
#         push!(trajectory, sum([create_ase_atoms(p, get_positions(config, i)) for i=1:size(positions[1])[3]]))
#     end

#     ase.io.write(file_name, trajectory)
# end

end # module