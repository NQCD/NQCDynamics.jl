"""
    QuantisedAtomic
This is a simplier version of "QuantisedDiatomic.jl", for dealing with an atom rather than a molecule.
This module exports two user facing functions:
- `generate_configurations`
    Creates a set of velocities and positions for an atom with specified
    translation energy
- `quantise_diatomic`
    Obtains translational energies for an atom with a given set of velocities and positions.
"""
module QuantisedAtomic

using LinearAlgebra: norm, normalize!, nullspace, cross, I
using UnitfulAtomic: austrip

using NQCDynamics: Simulation, Calculators, DynamicsUtils
using NQCBase: Atoms, PeriodicCell, InfiniteCell
using NQCModels: NQCModels
using NQCModels.AdiabaticModels: AdiabaticModel

export generate_configurations
#export quantise_diatomic

struct SurfaceParameters{T}
    reduced_mass::T
    atom_index
    slab::Matrix{T}
    height::T
    surface_normal::Vector{T}
    orthogonal_vectors::Matrix{T}
end

function SurfaceParameters(masses::AbstractVector, atom_index, slab::Matrix, height, surface_normal)
    μ = masses[atom_index]
    orthogonal_vectors = nullspace(reshape(surface_normal, 1, :)) # Plane parallel to surface
    if size(slab, 2) == 0
        surface_top = 0
    else
        surface_top = maximum(slab[3,:])
    end
    return SurfaceParameters{eltype(masses)}(μ, atom_index, slab, height + surface_top,
        surface_normal, orthogonal_vectors)
end

struct GenerationParameters{T}
    direction::Vector{T}
    translational_energy::T
    samples::Int
end

Base.broadcastable(p::SurfaceParameters) = Ref(p)
Base.broadcastable(p::GenerationParameters) = Ref(p)

"""
    generate_configurations(sim, ν, J; samples=1000, height=10, normal_vector=[0, 0, 1],
        translational_energy=0, direction=[0, 0, -1], position=[0, 0, height])
Generate positions and momenta for given quantum numbers
`translational_energy`, `direction` and `position` specify the kinetic energy in
a specific direction with the molecule placed with centre of mass at `position`.
Keyword arguments `height` and `normal_vector` become relevant if the potential
requires specific placement of the molecule.
These allow the molecule to be placed at a distance `height` in the direction
`normal_vector` when performing potential evaluations.
"""
function generate_configurations(sim;
    samples=1000,
    height=10.0,
    surface_normal=[0, 0, 1.0],
    translational_energy=0.0,
    direction=[0, 0, -1.0],
    atom_index=1,
    r=zeros(size(sim)))

    _, slab = separate_slab_and_molecule(atom_index, r)
    surface = SurfaceParameters(sim.atoms.masses, atom_index, slab, austrip(height), surface_normal)
    generation = GenerationParameters(direction, austrip(translational_energy), samples)

    # hack - bond is a dummy here (as well as myself)
    # to ensure we get samples amount of configurations (not just one)
    bonds = zeros(samples)

    @info "Generating the requested configurations..."
    configure_atomic.(sim, bonds, surface, generation)
end

function configure_atomic(sim,bond, surface::SurfaceParameters, generation::GenerationParameters)
    r = zeros(3,1)
    v = zeros(3,1)

    r[1,1] = bond
 
    position_above_surface!(r, surface.height, sim.cell)
    apply_translational_impulse!(v, sim.atoms.masses[surface.atom_index], generation.translational_energy, generation.direction)

    v, r
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

function position_above_surface!(r, height, cell::PeriodicCell)
    r[3,:] .+= height
    a1 = cell.vectors[:,1]
    a2 = cell.vectors[:,2]
    displacement = rand()*a1+rand()*a2
    r .+= displacement
    return nothing
end

function position_above_surface!(r, height, ::InfiniteCell)
    r[3,:] .+= height
    return nothing
end

function apply_translational_impulse!(v, masses, translational_energy, direction)
    velocity_impulse = velocity_from_energy(masses, translational_energy)
    v .+= velocity_impulse .* normalize!(direction)
end

function velocity_from_energy(masses, energy)
    m = sum(masses)
    sqrt(2austrip(energy) / m)
end

reduced_mass(atoms::Atoms) = reduced_mass(atoms.masses)
reduced_mass(m::AbstractVector) = m[1]*m[2]/(m[1]+m[2])

end # module