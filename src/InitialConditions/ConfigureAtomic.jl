"""
    ConfigureAtomic
This is a simplier version of "QuantisedDiatomic.jl", for dealing with an atom rather than a molecule.
This module exports one user facing function:
- `generate_configurations`
    Creates a set of velocities and positions for an atom with specified
    translation energy
"""
module ConfigureAtomic

using LinearAlgebra: norm, normalize!, nullspace, cross, I, dot
using UnitfulAtomic: austrip

using NQCDynamics: Simulation, Calculators, DynamicsUtils
using NQCBase: Atoms, PeriodicCell, InfiniteCell
using NQCModels: NQCModels
using NQCModels.AdiabaticModels: AdiabaticModel
using NQCDynamics.InitialConditions: QuantisedDiatomic
using Distributions: Uniform
export generate_configurations


"""
    generate_configurations(sim, ν, J; samples=1000, height=10, normal_vector=[0, 0, 1],
        translational_energy=0, direction=[0, 0, -1], position=[0, 0, height])
Generate positions and momenta for given quantum numbers
`translational_energy`, `direction` and `position` specify the kinetic energy in
a specific direction with the atom at `position`.
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
    incidence_angle=0.0, # Degrees
    azimuthal_angle=0.0,
    atom_index=[1],
    r=zeros(size(sim)))

    θ_i = zeros(samples)
    θ_i .= incidence_angle


    θ_j = zeros(samples)
    θ_j .= azimuthal_angle


    direction = set_incidence_angle.(θ_i,θ_j)

    _, slab = QuantisedDiatomic.separate_slab_and_molecule(atom_index, r)
    environment = QuantisedDiatomic.EvaluationEnvironment(atom_index, size(sim), slab, austrip(height), surface_normal)
    #surface = SurfaceParameters(sim.atoms.masses, atom_index, slab, austrip(height), surface_normal)
    generation = QuantisedDiatomic.GenerationParameters.(direction, austrip(translational_energy), samples)

    # hack - bond is a dummy here (as well as myself)
    # to ensure we get samples amount of configurations (not just one)
    bonds = zeros(samples)

    @info "Generating the requested configurations..."
    QuantisedDiatomic.configure_atomic.(sim, bonds, environment, generation)
    
end

"""
        set_incidence_direction(incidence_angle)
        Converts inputted angle in degrees to incidence direction vector for scattering simulations
"""
function set_incidence_angle(incidence_angle, azimuthal_angle)

    # Rotation in z to get incidence angle
    θ_i = deg2rad(incidence_angle)
    rot_mat = [
    1 0 0 
    0 cos(θ_i) -sin(θ_i)
    0 sin(θ_i) cos(θ_i)
    ]

    inc_dir = rot_mat * [0, 0, -1]

    # Random rotation in z
    #θ_2 = rand(Uniform(0.,2π),1,1)

    θ_2 = deg2rad(azimuthal_angle) 

    rot2_mat = [
    cos(θ_2) -sin(θ_2) 0
    sin(θ_2) cos(θ_2) 0
    0 0 1 
    ]


    
    inc_dir = rot2_mat * inc_dir

    inc_dir
end


function angle_to_surface_normal(sim::Simulation,v::Matrix;
    surface_normal=[0, 0, 1],
    atom_index = [1]
    )


    v_atom = v[:,atom_index][:]



    θ = atand(norm(cross(v_atom,surface_normal)),dot(v_atom,surface_normal))

    θ
end

reduced_mass(atoms::Atoms) = reduced_mass(atoms.masses)
reduced_mass(m::AbstractVector) = m[1]*m[2]/(m[1]+m[2])

end # module