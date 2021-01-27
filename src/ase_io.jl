
using .PyCall
using Unitful
using UnitfulAtomic

export read_system
export write_trajectory
export extract_parameters_and_positions
export create_ase_atoms
export read_trajectory

function read_system(file::String)
    io = pyimport("ase.io")
    atoms = io.read(file)

    extract_parameters_and_positions(atoms)
end

function extract_parameters_and_positions(ase_atoms::PyObject)
    Cell(ase_atoms), Atoms(ase_atoms), positions(ase_atoms)
end

Atoms(ase_atoms::PyObject) = Atoms{Float64}(Symbol.(ase_atoms.get_chemical_symbols()))

positions(ase_atoms::PyObject) = austrip.(ase_atoms.get_positions()'u"Å")

function Cell(ase_atoms::PyObject)
    if all(ase_atoms.cell.array .== 0)
        return InfiniteCell()
    else
        return PeriodicCell{Float64}(austrip.(ase_atoms.cell.array'u"Å"), ase_atoms.pbc)
    end
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

function read_trajectory(file_name::String)
    io = pyimport("ase.io")
    start = io.read(file_name, 0)
    cell = Cell(start)
    atoms = Atoms(start)

    trajectory = io.iread(file_name)
    R = []
    for ase_atoms in trajectory
        push!(R, positions(ase_atoms))
    end
    (cell, atoms, R)
end
