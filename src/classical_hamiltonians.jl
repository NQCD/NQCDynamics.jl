
evaluate_hamiltonian(sim::AbstractSimulation, u) =
    evaluate_hamiltonian(sim, DynamicsUtils.get_velocities(u), DynamicsUtils.get_positions(u))

function evaluate_hamiltonian(sim::AbstractSimulation, v, r)
    k = evaluate_kinetic_energy(sim.atoms.masses, v)
    e = evaluate_potential_energy(sim, r)
    k + e
end

evaluate_kinetic_energy(masses, v) = sum(masses' .* v.^2)/2
function evaluate_kinetic_energy(masses, v::RingPolymerArray)
    E = 0.0
    @views for i in axes(v, 3)
        E += evaluate_kinetic_energy(masses, v[:,:,i])
    end
    E
end

function evaluate_potential_energy(sim::Simulation, R)
    Calculators.evaluate_potential!(sim.calculator, R)
    sim.calculator.potential[1]
end

function evaluate_potential_energy(sim::RingPolymerSimulation, R)
    Calculators.evaluate_potential!(sim.calculator, R)
    sum(sim.calculator.potential)[1] + get_spring_energy(sim, R)
end

function evaluate_potential_energy(sim::RingPolymerSimulation, R::RingPolymerArray)
    if R.normal
        transform!(sim.beads, R)
        potential = evaluate_potential_energy(sim, R.data)
        transform!(sim.beads, R)
    else
        potential = evaluate_potential_energy(sim, R.data)
    end
    return potential
end

"""
    get_spring_energy(system::RingPolymerSimulation, R)
    
Calculate the ring polymer spring potential.
"""
function get_spring_energy(sim::RingPolymerSimulation, R)
    E = 0.0
    for bead=1:length(sim.beads)-1
        for i in sim.beads.quantum_atoms # Only for quantum nuclei
            for j=1:sim.DoFs
                E += sim.atoms.masses[i] * (R[j, i, bead] - R[j, i, bead+1])^2
            end
        end
    end
    for i in sim.beads.quantum_atoms
        for j=1:sim.DoFs
            E += sim.atoms.masses[i] * (R[j, i, length(sim.beads)] - R[j, i, 1])^2
        end
    end
    E * sim.beads.ω_n^2/2
end
