"""
    QuantisedDiatomic

This module exports two user facing functions:

- `generate_configurations`
    Creates a set of velocities and positions for diatomic molecule with specified
    vibrational `ν` and rotational `J` quantum numbers.

- `quantise_diatomic`
    Obtains vibrational `ν` and rotational `J` quantum numbers for a diatomic molecule
    with a given set of velocities and positions.

The central concept of this module is the EBK procedure which is nicely detailed here:
[Larkoski2006](@cite)

Inspired by VENUS96: [Hase1996](@cite)
"""
module QuantisedDiatomic

using QuadGK: QuadGK
using Roots: Roots
using Rotations: Rotations
using LsqFit: LsqFit
using LinearAlgebra: norm, normalize!, nullspace, cross, I
using UnitfulAtomic: austrip
using Distributions: Uniform

using NonadiabaticMolecularDynamics: Simulation, Calculators, DynamicsUtils
using NonadiabaticDynamicsBase: Atoms, PeriodicCell, InfiniteCell

export generate_configurations
export quantise_diatomic

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
function generate_configurations(sim, ν, J; samples=1000, height=10, normal_vector=[0, 0, 1],
    translational_energy=0, direction=[0, 0, -1])
    E, bounds = extract_energy_and_bounds(sim, ν, J; height=height, normal_vector=normal_vector)
    μ = reduced_mass(sim.atoms)

    rotational(r) = J*(J+1) / (2μ*r^2)
    Vᵣ(r) = rotational(r) + calculate_diatomic_energy(sim, r; height=height, normal_vector=normal_vector)
    radial_momentum(r) = sqrt(2μ * (E - Vᵣ(r)))

    bonds, momenta = select_random_bond_lengths(E, bounds, radial_momentum, μ, samples)

    configure_diatomic.(sim, bonds, momenta, J, μ, Ref(direction), translational_energy,height)
end

"""
    extract_energy_and_bounds(sim, ν, J; height=10, normal_vector=[0, 0, 1])

Returns the energy and turning points associated with the specified quantum numbers
"""
function extract_energy_and_bounds(sim, ν, J; height=10, normal_vector=[0, 0, 1])

    k = calculate_force_constant(sim; height=height, normal_vector=normal_vector)
    μ = reduced_mass(sim.atoms)
    ω = sqrt(k / μ)
    E_guess = (ν + 1/2) * ω

    "Target function to optimise. Find energy `E` where vibrational numbers `ν` match."
    function target(E)
        ν_tmp, _, _ = extract_quantum_numbers(sim, μ, E, J; height=height, normal_vector=normal_vector)
        ν_tmp - ν
    end

    E = Roots.find_zero(target, E_guess)
    _, _, bounds = extract_quantum_numbers(sim, μ, E, J; height=height, normal_vector=normal_vector)

    E, bounds 
end

"""
    quantise_diatomic(sim::Simulation, v::Matrix, r::Matrix;
        height=10, normal_vector=[0, 0, 1])

Quantise the vibrational and rotational degrees of freedom for the specified
positions and velocities

When evaluating the potential, the molecule is moved to `height` in direction `normal_vector`.
If the potential is independent of centre of mass position, this has no effect.
Otherwise, be sure to modify these parameters to give the intended behaviour.
"""
function quantise_diatomic(sim::Simulation, v::Matrix, r::Matrix;
    height=10, normal_vector=[0, 0, 1])
    # Select only the first two atoms
    if length(sim.atoms) > 2
        sim = Simulation(sim.atoms[1:2], sim.calculator.model; cell=sim.cell)
        v = v[:,1:2]
        r = r[:,1:2]
    elseif length(sim.atoms) == 1
        throw(ArgumentError("Only a single atom in the simulation?"))
    end

    r_com = subtract_centre_of_mass(r, sim.atoms.masses)
    v_com = subtract_centre_of_mass(v, sim.atoms.masses)
    p_com = v_com .* sim.atoms.masses'

    k = DynamicsUtils.classical_kinetic_energy(sim, v_com)
    p = calculate_diatomic_energy(sim, bond_length(r_com); height=height, normal_vector=normal_vector)
    E = k + p

    L = total_angular_momentum(r_com, p_com)
    J = (sqrt(1+4*L^2) - 1) / 2 # L^2 = J(J+1)ħ^2

    μ = reduced_mass(sim.atoms)
    ν, J, bounds = extract_quantum_numbers(sim, μ, E, J; height=height, normal_vector=normal_vector)
    ν, J
end

"""
Calculate the vibrational and rotational quantum numbers using EBK
"""
function extract_quantum_numbers(sim::Simulation, μ::Real, E::Real, J::Real; height=10, normal_vector=[0, 0, 1])

    rotational(r) = J*(J+1) / (2μ*r^2)
    Vᵣ(r) = rotational(r) + calculate_diatomic_energy(sim, r)
    radial_momentum(r) = sqrt(2μ * (E - Vᵣ(r)))

    bounds = Roots.find_zeros(r -> E - Vᵣ(r), 0.0, 10.0)
    if length(bounds) != 2
        throw(ArgumentError("Unable to determine bounds for EBK quantisation."
        * " Perhaps the chosen quantum numbers are unreasonable?
        Bounds obtained = $bounds"))
    end
    integral, err = QuadGK.quadgk(radial_momentum, bounds...)
    ν = integral / π - 1/2

    ν, J, bounds
end

"""
Pick a random bond length and corresponding radial momentum that matches the
radial probability distribution.
"""
function select_random_bond_lengths(E, bounds, probability_function, μ, samples)
    P0 = sqrt(1e-4 * 2μ*E)
    distribution = Uniform(bounds...)
    bonds = zeros(samples)
    momenta = zeros(samples)
    for i=1:samples
        keep_going = true
        while keep_going
            r, P = pick_radius(distribution, probability_function, 1)
            if rand() < P0 / P
                bonds[i] = r
                momenta[i] = rand([-1, 1]) * P
                keep_going = false
            end
        end
    end
    bonds, momenta
end

function pick_radius(distribution, func, count)
    r = rand(distribution)
    try
        P = func(r)
        return r, P
    catch e
        if e isa DomainError
            if count > 100
                throw(ArgumentError("Something is seriously wrong, a rejected bond length was generated 100 times in a row?
                    Check the model."))
            else
                return pick_radius(distribution, func, count+1)
            end
        else
            throw(e)
        end
    end
end

"""
Randomly orient molecule in space for a given bond length and radial momentum
"""
function configure_diatomic(sim, bond, momentum, J, μ, direction, translational_energy,height)
    r = zeros(3,2)
    v = zeros(3,2)
    r[1,1] = bond
    v[1,1] = momentum ./ μ
    r = subtract_centre_of_mass(r, sim.atoms.masses)
    v = subtract_centre_of_mass(v, sim.atoms.masses)

    L = sqrt(J * (J + 1))

    random = 2π*rand()
    ω = L .* [0, sin(random), cos(random)] ./ (μ * bond^2)

    v[:,1] .+= cross(ω, r[:,1])
    v[:,2] .+= cross(ω, r[:,2])

    apply_random_rotation!(v, r)
    position_above_surface!(r, height, sim.cell)
    apply_translational_impulse!(v, sim.atoms.masses, translational_energy, direction)

    v, r
end

"""
Randomly rotate each column of two 3*N matrix, same rotation for all columns.
"""
function apply_random_rotation!(x, y)
    rotation = rand(Rotations.RotMatrix{3})
    for (col1, col2) in zip(eachcol(x), eachcol(y))
        col1 .= rotation * col1
        col2 .= rotation * col2
    end
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

function velocity_from_energy(masses, energy)
    m = sum(masses)
    sqrt(2austrip(energy) / m)
end

function apply_translational_impulse!(v, masses, translational_energy, direction)
    velocity_impulse = velocity_from_energy(masses, translational_energy)
    v .+= velocity_impulse .* normalize!(direction)
end

"""
    calculate_diatomic_energy(sim::Simulation, bond_length::Real;
        height=10, normal_vector=[0, 0, 1])

Returns potential energy of diatomic with `bond_length` at `height` from surface.

Orients molecule parallel to the surface at the specified height, the surface is assumed
to intersect the origin.
This requires that the model implicitly provides the surface, or works fine without one.
"""
function calculate_diatomic_energy(sim::Simulation, bond_length::Real;
    height=10, normal_vector=[0, 0, 1])

    orthogonals = nullspace(reshape(normal_vector, 1, :)) # Plane parallel to surface
    R = normal_vector .* height .+ orthogonals .* bond_length ./ sqrt(2)

    Calculators.evaluate_potential!(sim.calculator, R)

    sim.calculator.potential[1]
end

"""
    calculate_force_constant(sim::Simulation;
        height=10.0, bond_length=2.5, normal_vector=[0.0, 0.0, 1.0])

Evaluate energy for different bond lengths and identify force constant. 

This uses LsqFit.jl to fit a harmonic potential with optimised minimum.
"""
function calculate_force_constant(sim::Simulation;
    height=10.0, bond_length=2.5, normal_vector=[0.0, 0.0, 1.0])

    harmonic(f, x, p) = (@. f = p[1]*(x - p[2])^2/2)
    function jacobian(J, x, p)
        @. J[:,1] = (x - p[2])^2/2
        @. J[:,2] = -p[1] * (x - p[2])
    end

    bond_lengths = 0.5:0.01:4.0
    V = calculate_diatomic_energy.(sim, bond_lengths; height=height, normal_vector=normal_vector)

    fit = LsqFit.curve_fit(harmonic, bond_lengths, V, [1.0, bond_length]; lower=[0.0, 0.0], inplace=true)
    fit.param[1]
end

subtract_centre_of_mass(x, m) = x .- centre_of_mass(x, m)
@views centre_of_mass(x, m) = (x[:,1].*m[1] .+ x[:,2].*m[2]) ./ (m[1]+m[2])
reduced_mass(atoms::Atoms) = reduced_mass(atoms.masses)
reduced_mass(m::AbstractVector) = m[1]*m[2]/(m[1]+m[2])
@views bond_length(r) = norm(r[:,1] .- r[:,2])
@views total_angular_momentum(r, p) = norm(cross(r[:,1], p[:,1]) + cross(r[:,2], p[:,2]))
@views total_moment_of_inertia(r, m) = norm(r[:,1])*m[1] + norm(r[:,2])*m[2]

function moment_of_inertia_tensor(r, m)
    tensor = zeros(3,3)
    for atom in axes(r, 2)
        absr = norm(r[:,atom])^2
        tensor += m[atom] * absr * I
        for i in axes(r, 1)
            for j in axes(r, 1)
                tensor[i,j] -= m[atom] * r[i,atom]*r[j,atom]
            end
        end
    end
    tensor
end

end # module
