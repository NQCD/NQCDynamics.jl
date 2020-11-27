module IO

using PyCall
using Unitful
using ..Systems

export read_system

function read_system(file::String)
    io = pyimport("ase.io")
    atoms = io.read(file)

    cell = Cell(atoms.cell.data, u"Å")
    atom_types = Symbol.(atoms.get_chemical_symbols())
    p = SystemParameters(cell, atom_types)

    R = vcat(atoms.get_positions()'...)
    P = zero(R)
    z = Phasespace(R, P)
    
    System(p, z)
end

end # module