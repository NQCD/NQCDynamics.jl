import AtomsBase
using StaticArraysCore: SVector

function NQCBase.Atoms(system::AtomsBase.AbstractSystem)
    return NQCBase.Atoms(Symbol.(AtomsBase.atomic_symbol(system, :))) # AtomsBase requires an index to be given and doesn't return as Symbol any more.
end

function Cell(system::AtomsBase.AbstractSystem)
    if isa(AtomsBase.cell(system), AtomsBase.IsolatedCell)
        return NQCBase.InfiniteCell()
    else
        box = AtomsBase.cell_vectors(system)
        cell = PeriodicCell(reduce(hcat, box))
        NQCBase.set_periodicity!(cell, vcat(AtomsBase.periodicity(system)...))
        return cell
    end
end

function Position(system::AtomsBase.AbstractSystem)
    r = AtomsBase.position(system, :)
    output = zeros(AtomsBase.n_dimensions(system), Base.length(system))
    for i in axes(output, 2)
        for j in axes(output, 1)
            output[j, i] = austrip(r[i][j])
        end
    end
    return output
end

function Velocity(system::AtomsBase.AbstractSystem)
    v = AtomsBase.velocity(system, :)
    output = zeros(AtomsBase.n_dimensions(system), Base.length(system))
    for i in axes(output, 2)
        for j in axes(output, 1)
            output[j, i] = austrip(v[i][j])
        end
    end
    return output
end

function AtomsBase.cell_vectors(cell::PeriodicCell)
    S = size(cell.vectors, 2)
    return SVector{S}(auconvert.(u"Å", vec) for vec in eachcol(cell.vectors))
end


function System(atoms::NQCBase.Atoms, position::AbstractMatrix, cell::AbstractCell=InfiniteCell())
    output_atoms = AtomsBaseAtoms(atoms, position)
    return build_system(output_atoms, cell)
end
function System(structure::NQCBase.Structure)
    return System(structure.atoms, structure.positions, structure.cell)
end

function System(atoms::NQCBase.Atoms, position::AbstractMatrix, velocity::AbstractMatrix, cell::AbstractCell=InfiniteCell())
    output_atoms = AtomsBaseAtoms(atoms, position, velocity)
    return build_system(output_atoms, cell)
end

function AtomsBaseAtoms(atoms::NQCBase.Atoms, position::AbstractMatrix)
    if length(atoms) != size(position, 2)
        @error atoms position
        error("The provided `Atoms` do not match the `position` array.")
    end

    output_atoms = AtomsBase.Atom[]
    sizehint!(output_atoms, length(atoms))
    for i in axes(position,2)
        r = auconvert.(u"Å", position[:,i])
        push!(output_atoms, AtomsBase.Atom(atoms.numbers[i], r))
    end
    return output_atoms
end

function AtomsBaseAtoms(atoms::NQCBase.Atoms, position::AbstractMatrix, velocity::AbstractMatrix)
    if length(atoms) != size(position, 2)
        @error atoms position
        error("The provided `Atoms` do not match the `position` array.")
    end
    if length(atoms) != size(velocity, 2)
        @error atoms velocity
        error("The provided `Atoms` do not match the `velocity` array.")
    end

    output_atoms = AtomsBase.Atom[]
    sizehint!(output_atoms, length(atoms))
    for i in axes(position,2)
        r = auconvert.(u"Å", position[:,i])
        v = auconvert.(u"Å/ps", velocity[:,i])
        push!(output_atoms, AtomsBase.Atom(atoms.numbers[i], r, v))
    end
    return output_atoms
end

function build_system(atoms, cell::InfiniteCell)
    return AtomsBase.isolated_system(atoms)
end

function build_system(atoms, cell::PeriodicCell)
    box = AtomsBase.cell_vectors(cell)
    bc = (cell.periodicity...,) # PBC are now a Tuple{Bool} to make everything harder to use.
    return AtomsBase.atomic_system(atoms, box, bc)
end

function Trajectory(
    atoms::NQCBase.Atoms,
    position::Vector{<:AbstractMatrix},
    velocity::Vector{<:AbstractMatrix},
    cell::AbstractCell=InfiniteCell()
    )

    trajectory = AtomsBase.FlexibleSystem[]
    sizehint!(trajectory, length(position))

    for i in eachindex(position, velocity)
        push!(trajectory, System(atoms, position[i], velocity[i], cell))
    end

    return trajectory
end
