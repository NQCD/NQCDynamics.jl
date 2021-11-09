using LinearAlgebra: eigvecs, diag, diagind, dot
using StatsBase: sample, Weights

using NonadiabaticMolecularDynamics.Calculators: RingPolymerDiabaticCalculator
using NonadiabaticMolecularDynamics: RingPolymerSimulation, RingPolymers

function RingPolymerSimulation{FSSH}(atoms::Atoms{S,T}, model::Model, n_beads::Integer; rescaling=:standard, kwargs...) where {S,T}
    RingPolymerSimulation(atoms, model, FSSH{T}(NonadiabaticModels.nstates(model), rescaling), n_beads; kwargs...)
end

function DynamicsMethods.motion!(du, u, sim::RingPolymerSimulation{<:SurfaceHopping}, t)
    dr = DynamicsUtils.get_positions(du)
    dv = DynamicsUtils.get_velocities(du)
    dσ = DynamicsUtils.get_quantum_subsystem(du)

    r = DynamicsUtils.get_positions(u)
    v = DynamicsUtils.get_velocities(u)
    σ = DynamicsUtils.get_quantum_subsystem(u)

    set_state!(u, sim.method.state) # Make sure the state variables match, 

    DynamicsUtils.velocity!(dr, v, r, sim, t)
    Calculators.update_electronics!(sim.calculator, r)
    acceleration!(dv, v, r, sim, t, sim.method.state)
    DynamicsUtils.apply_interbead_coupling!(dv, r, sim)
    DynamicsUtils.set_quantum_derivative!(dσ, v, σ, sim)
end

function get_hopping_nonadiabatic_coupling(sim::RingPolymerSimulation{<:FSSH})
    sim.calculator.centroid_nonadiabatic_coupling
end

function get_hopping_velocity(::RingPolymerSimulation{<:FSSH}, u)
    RingPolymers.get_centroid(DynamicsUtils.get_velocities(u))
end

function get_hopping_eigenvalues(sim::RingPolymerSimulation{<:FSSH})
    sim.calculator.centroid_eigenvalues
end

function calculate_potential_energy_change(calc::RingPolymerDiabaticCalculator, new_state::Integer, current_state::Integer)
    return calc.centroid_eigenvalues[new_state] - calc.centroid_eigenvalues[current_state]
end

function perform_rescaling!(sim::RingPolymerSimulation{<:FSSH}, velocity, velocity_rescale, new_state, old_state)
    d = get_hopping_nonadiabatic_coupling(sim)
    for i in range(sim.atoms)
        coupling = [d[j,i][new_state, old_state] for j=1:ndofs(sim)]
        velocity[:,i,:] .-= velocity_rescale[i] .* coupling ./ sim.atoms.masses[i]
    end
    return nothing
end

function DynamicsUtils.classical_hamiltonian(sim::RingPolymerSimulation{<:FSSH}, u)
    kinetic = DynamicsUtils.classical_kinetic_energy(sim, DynamicsUtils.get_velocities(u))
    spring = RingPolymers.get_spring_energy(sim.beads, sim.atoms.masses, DynamicsUtils.get_positions(u))

    Calculators.evaluate_potential!(sim.calculator, DynamicsUtils.get_positions(u))
    Calculators.eigen!(sim.calculator)
    potential = sum(bead[u.state] for bead in sim.calculator.eigenvalues)
    return kinetic + spring + potential
end
