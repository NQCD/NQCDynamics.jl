using QuadGK
using Distributions
using Roots

"""
    Calculate the vibrational and rotational quantum numbers using EBK
"""
function extract_quantum_numbers(sim::Simulation, r::Matrix, p::Matrix)
    @assert length(sim.atoms) == 2
    r_com = subtract_centre_of_mass(r, sim.atoms.masses)
    p_com = subtract_centre_of_mass(p, sim.atoms.masses)

    μ = reduced_mass(sim.atoms)

    v_com = p_com ./ sim.atoms.masses'
    E = evaluate_hamiltonian(sim, v_com, r)

    L = total_angular_momentum(r_com, p_com)

    extract_quantum_numbers(sim, μ, E, L)
end

function extract_quantum_numbers(sim::Simulation, μ::Real, E::Real, L::Real)

    J = (sqrt(1+4*L^2) - 1) / 2 # L^2 = J(J+1)ħ^2

    rotational(r) = J*(J+1) / (2μ*r^2)
    Vᵣ(r) = rotational(r) + calculate_diatomic_energy(sim, r)

    radial_momentum(r) = sqrt(2μ * (E - Vᵣ(r)))

    bounds = find_zeros(r -> E - Vᵣ(r), 0.0, 50.0)
    integral, err = quadgk(radial_momentum, bounds...)
    ν = integral / π - 1/2

    ν, J, bounds, integral
end

"""
    Identify the energy and turning points associated with the specified quantum numbers
"""
function extract_energy_and_bounds(sim, ν, J)
    L = sqrt(J * (J + 1))

    k = calculate_force_constant(sim)
    μ = reduced_mass(sim.atoms)
    ω = sqrt(k / μ)
    E_guess = (ν + 1/2) * ω

    function target(E)
        ν_tmp, J_tmp, bounds, integral = extract_quantum_numbers(sim, μ, E, L)
        ν_tmp - ν
    end

    E = find_zero(target, E_guess)
    ν_tmp, J_tmp, bounds, integral = extract_quantum_numbers(sim, μ, E, L)

    E, bounds, integral
end

"""
    Generate positions and momenta for given quantum numbers
"""
function generate_configuration(sim, ν, J; samples=1000)
    E, bounds, integral = extract_energy_and_bounds(sim, ν, J)
    μ = reduced_mass(sim.atoms)

    rotational(r) = J*(J+1) / (2μ*r^2)
    Vᵣ(r) = rotational(r) + calculate_diatomic_energy(sim, r)
    radial_momentum(r) = sqrt(2μ * (E - Vᵣ(r)))

    bonds = select_random_bond_lengths(E, bounds, radial_momentum, μ, samples)

end

function select_random_bond_lengths(E, bounds, probability_function, μ, samples)
    P0 = sqrt(1e-4 * 2μ*E)
    distribution = Uniform(bounds...)
    bonds = zeros(samples)
    for i=1:samples
        keep_going = true
        while keep_going
            r = rand(distribution)
            P = probability_function(r)
            if rand() < P0 / P
                bonds[i] = r
                keep_going = false
            end
        end
    end
    bonds
end

end