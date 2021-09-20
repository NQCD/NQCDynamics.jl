
evaluate_hamiltonian(sim::AbstractSimulation, u) =
    evaluate_hamiltonian(sim, DynamicsUtils.get_velocities(u), DynamicsUtils.get_positions(u))

function evaluate_hamiltonian(sim::AbstractSimulation, v, r)
    k = evaluate_kinetic_energy(sim.atoms.masses, v)
    e = evaluate_potential_energy(sim, r)
    k + e
end

evaluate_kinetic_energy(masses, v) = sum(masses' .* v.^2)/2
function evaluate_kinetic_energy(masses, v::RingPolymers.RingPolymerArray)
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
    sum(sim.calculator.potential)[1] + RingPolymers.get_spring_energy(sim.beads, sim.atoms.masses, R)
end

function evaluate_potential_energy(sim::RingPolymerSimulation, R::RingPolymers.RingPolymerArray)
    if R.normal
        RingPolymers.transform!(sim.beads, R)
        potential = evaluate_potential_energy(sim, R.data)
        RingPolymers.transform!(sim.beads, R)
    else
        potential = evaluate_potential_energy(sim, R.data)
    end
    return potential
end
