using NonadiabaticMolecularDynamics: RingPolymers

function RingPolymerSimulation{Ehrenfest}(atoms::Atoms{S,T}, model::Model, n_beads::Integer; kwargs...) where {S,T}
    RingPolymerSimulation(atoms, model, Ehrenfest{T}(NonadiabaticModels.nstates(model)), n_beads; kwargs...)
end

function DynamicsMethods.motion!(du, u, sim::RingPolymerSimulation{<:Ehrenfest}, t)
    dr = DynamicsUtils.get_positions(du)
    dv = DynamicsUtils.get_velocities(du)
    dσ = DynamicsUtils.get_quantum_subsystem(du)

    r = DynamicsUtils.get_positions(u)
    v = DynamicsUtils.get_velocities(u)
    σ = DynamicsUtils.get_quantum_subsystem(u)

    DynamicsUtils.velocity!(dr, v, r, sim, t)
    Calculators.update_electronics!(sim.calculator, r)
    acceleration!(dv, v, r, sim, t, σ)
    DynamicsUtils.apply_interbead_coupling!(dv, r, sim)
    DynamicsUtils.set_quantum_derivative!(dσ, v, σ, sim)
end

function Estimators.diabatic_population(sim::RingPolymerSimulation{<:Ehrenfest}, u)
    Calculators.evaluate_centroid_potential!(sim.calculator, DynamicsUtils.get_positions(u))
    Calculators.centroid_eigen!(sim.calculator)
    U = sim.calculator.centroid_eigenvectors
    σ = DynamicsUtils.get_quantum_subsystem(u)
    return real.(diag(U * σ * U'))
end

function DynamicsUtils.classical_hamiltonian(sim::RingPolymerSimulation{<:Ehrenfest}, u)
    v = DynamicsUtils.get_velocities(u)
    r = DynamicsUtils.get_positions(u)
    kinetic = DynamicsUtils.classical_kinetic_energy(sim, v)
    spring = RingPolymers.get_spring_energy(sim.beads, sim.atoms.masses, r)

    Calculators.evaluate_potential!(sim.calculator, r)
    Calculators.eigen!(sim.calculator)
    population = Estimators.adiabatic_population(sim, u)
    potential = sum([dot(population, eigs) for eigs in sim.calculator.eigenvalues])

    return kinetic + potential + spring
end
