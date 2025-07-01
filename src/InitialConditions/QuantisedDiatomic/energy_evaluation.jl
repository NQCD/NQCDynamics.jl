
struct EvaluationEnvironment{T}
    molecule_indices::Vector{Int}
    system_size::Tuple{Int,Int}
    other_atoms::Matrix{T}
    offset::Vector{T}
    orthogonal_vectors::Matrix{T}
end

Base.broadcastable(p::EvaluationEnvironment) = Ref(p)

function EvaluationEnvironment(molecule_indices, system_size, other_atoms, height, surface_normal)
    normalize!(surface_normal) # Normalise surface normal vector
    # If surface atoms are explicitly included, find the height of the surface along `surface_normal` and adjust `height` accordingly.  
    if size(other_atoms, 2) != 0
        height += reduce(max, mapslices(norm, other_atoms, dims=1))
    end
    offset = surface_normal .* height
    orthogonal_vectors = nullspace(reshape(surface_normal, 1, :)) # Plane parallel to surface
    return EvaluationEnvironment(molecule_indices, system_size, other_atoms, offset, orthogonal_vectors)
end

"""
    separate_slab_and_molecule(atom_indices, r)

Get the coordinates of the molecule and slab separately.
"""
function separate_slab_and_molecule(atom_indices, r)
    molecule = r[:,atom_indices]
    slab_indices = [i for i in axes(r, 2) if i ∉ atom_indices]
    slab = r[:,slab_indices]
    return (molecule, slab)
end

"""
    combine_slab_and_molecule(atom_indices, molecule, slab)

Revert the transformation `separate_slab_and_molecule`
"""
function combine_slab_and_molecule(atom_indices, molecule, slab)
    r = zeros(size(molecule, 1), size(molecule, 2) + size(slab, 2))
    slab_indices = [i for i in axes(r, 2) if i ∉ atom_indices]
    r[:,atom_indices] .= molecule
    r[:,slab_indices] .= slab
    return r
end

"""
    calculate_diatomic_energy(model::ClassicalModel, bond_length::Real;
        height=10, normal_vector=[0, 0, 1])

Returns potential energy of diatomic with `bond_length` at `height` from surface.

Orients molecule parallel to the surface at the specified height within the simulation cell, 
assuming the height has already been adjusted to include that of the surface. 

(this is checked in the EvaluationEnvironment constructor)
"""
function calculate_diatomic_energy(
    bond_length::Real, model::ClassicalModel, environment::EvaluationEnvironment
)
    r = assemble_evaluation_geometry(bond_length, environment)
    return NQCModels.potential(model, r)
end

function assemble_evaluation_geometry(bond_length::Real, environment::EvaluationEnvironment)
    molecule = build_molecule(bond_length, environment)
    r = combine_slab_and_molecule(environment.molecule_indices, molecule, environment.other_atoms)
    return r
end

function build_molecule(bond_length::Real, environment::EvaluationEnvironment)
    if environment.system_size == (1, 1)
        r = [bond_length;;]
    else
        r = environment.offset .+ environment.orthogonal_vectors .* bond_length ./ sqrt(2)
    end
    return r
end