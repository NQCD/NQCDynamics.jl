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
Am. J. Phys., Vol. 74, No. 7, July 2006

Inspired by VENUS96:
Hase, WL et al. 1996, VENUS96: A general chemical dynamics computer program,
Quantum Chemical Program Exchange Bulletin, vol. 16(4), pp. 43-43.
"""
module QuantisedDiatomic

using QuadGK
using Distributions
using Roots
using Rotations
using LsqFit
using LinearAlgebra
using ....NonadiabaticMolecularDynamics

export generate_configurations
export quantise_diatomic

"""
    generate_configurations(sim, ν, J; samples=1000)

Generate positions and momenta for given quantum numbers
"""
function generate_configurations(sim, ν, J; samples=1000)
    E, bounds = extract_energy_and_bounds(sim, ν, J)
    μ = reduced_mass(sim.atoms)

    rotational(r) = J*(J+1) / (2μ*r^2)
    Vᵣ(r) = rotational(r) + calculate_diatomic_energy(sim, r)
    radial_momentum(r) = sqrt(2μ * (E - Vᵣ(r)))

    bonds, momenta = select_random_bond_lengths(E, bounds, radial_momentum, μ, samples)

    configure_diatomic.(sim, bonds, momenta, J, μ)
end

"""
Identify the energy and turning points associated with the specified quantum numbers
"""
function extract_energy_and_bounds(sim, ν, J)

    k = calculate_force_constant(sim)
    μ = reduced_mass(sim.atoms)
    ω = sqrt(k / μ)
    E_guess = (ν + 1/2) * ω

    function target(E)
        ν_tmp, J_tmp, bounds = extract_quantum_numbers(sim, μ, E, J)
        ν_tmp - ν
    end

    E = find_zero(target, E_guess)
    ν_tmp, J_tmp, bounds = extract_quantum_numbers(sim, μ, E, J)

    E, bounds 
end

"""
    quantise_diatomic(sim::Simulation, v::Matrix, r::Matrix)

Quantise the vibrational and rotational degrees of freedom for the specified
positions and velocities
"""
function quantise_diatomic(sim::Simulation, v::Matrix, r::Matrix)
    @assert length(sim.atoms) == 2
    p = v .* sim.atoms.masses'
    r_com = subtract_centre_of_mass(r, sim.atoms.masses)
    p_com = subtract_centre_of_mass(p, sim.atoms.masses)

    μ = reduced_mass(sim.atoms)

    v_com = p_com ./ sim.atoms.masses'
    E = evaluate_hamiltonian(sim, v_com, r)

    L = total_angular_momentum(r_com, p_com)
    J = (sqrt(1+4*L^2) - 1) / 2 # L^2 = J(J+1)ħ^2

    ν, J, bounds = extract_quantum_numbers(sim, μ, E, J)
    ν, J
end

"""
Calculate the vibrational and rotational quantum numbers using EBK
"""
function extract_quantum_numbers(sim::Simulation, μ::Real, E::Real, J::Real)

    rotational(r) = J*(J+1) / (2μ*r^2)
    Vᵣ(r) = rotational(r) + calculate_diatomic_energy(sim, r)
    radial_momentum(r) = sqrt(2μ * (E - Vᵣ(r)))

    bounds = find_zeros(r -> E - Vᵣ(r), 0.0, 10.0)
    integral, err = quadgk(radial_momentum, bounds...)
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
            r = rand(distribution)
            P = probability_function(r)
            if rand() < P0 / P
                bonds[i] = r
                momenta[i] = rand([-1, 1]) * P
                keep_going = false
            end
        end
    end
    bonds, momenta
end

"""
Randomly orient molecule in space for a given bond length and radial momentum
"""
function configure_diatomic(sim, bond, momentum, J, μ)
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

    v, r
end

"""
Randomly rotate each column of two 3*N matrix, same rotation for all columns.
"""
function apply_random_rotation!(x, y)
    rotation = rand(RotMatrix{3})
    for (col1, col2) in zip(eachcol(x), eachcol(y))
        col1 .= rotation * col1
        col2 .= rotation * col2
    end
end

"""
    calculate_diatomic_energy(sim::Simulation, bond_length::T; height::T=10.0,
                               normal_vector::Vector{T}=[0.0, 0.0, 1.0]) where {T}

Returns potential energy of diatomic with `bond_length` at `height` from surface.

Orients molecule parallel to the surface at the specified height, the surface is assumed
to intersect the origin.
This requires that the model implicitly provides the surface, or works fine without one.
"""
function calculate_diatomic_energy(sim::Simulation, bond_length::Real;
                                   height::Real=10.0,
                                   normal_vector::Vector=[0.0, 0.0, 1.0])

    orthogonals = nullspace(reshape(normal_vector, 1, :)) # Plane parallel to surface
    R = normal_vector .* height .+ orthogonals .* bond_length ./ sqrt(2)

    Calculators.evaluate_potential!(sim.calculator, R)

    sim.calculator.potential[1]
end

"""
Evaluate energy for different bond lengths and identify force constant. 

This uses LsqFit.jl to fit a harmonic potential with optimised minimum.
"""
function calculate_force_constant(sim::Simulation;
                                  height::Real=10.0,
                                  bond_length=2.5,
                                  normal_vector::Vector=[0.0, 0.0, 1.0])

    harmonic(f, x, p) = (@. f = p[1]*(x - p[2])^2/2)
    function jacobian(J, x, p)
        @. J[:,1] = (x - p[2])^2/2
        @. J[:,2] = -p[1] * (x - p[2])
    end

    bond_lengths = 0.5:0.01:4.0
    V = calculate_diatomic_energy.(sim, bond_lengths; height=height, normal_vector=normal_vector)

    fit = curve_fit(harmonic, bond_lengths, V, [1.0, bond_length]; lower=[0.0, 0.0], inplace=true)
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