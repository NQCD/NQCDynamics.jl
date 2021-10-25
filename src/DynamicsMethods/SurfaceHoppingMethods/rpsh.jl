using LinearAlgebra: eigvecs, diag, diagind, dot
using StatsBase: sample, Weights

using NonadiabaticMolecularDynamics.Calculators: RingPolymerDiabaticCalculator
using NonadiabaticMolecularDynamics: RingPolymerSimulation, RingPolymers

function RingPolymerSimulation{FSSH}(atoms::Atoms{S,T}, model::Model, n_beads::Integer; kwargs...) where {S,T}
    RingPolymerSimulation(atoms, model, FSSH{T}(NonadiabaticModels.nstates(model)), n_beads; kwargs...)
end

function DynamicsMethods.DynamicsVariables(sim::RingPolymerSimulation{<:SurfaceHopping}, v, r, electronic::NonadiabaticDistributions.SingleState)
    n_states = NonadiabaticModels.nstates(sim.calculator.model)
    state = electronic.state
    if electronic.statetype === NonadiabaticDistributions.Diabatic()
        Calculators.evaluate_centroid_potential!(sim.calculator, r)
        U = eigvecs(sim.calculator.centroid_potential)

        diabatic_density = zeros(n_states, n_states)
        diabatic_density[state, state] = 1
        σ = U' * diabatic_density * U
        state = sample(Weights(diag(real.(σ))))
    else
        σ = zeros(n_states, n_states)
        σ[state, state] = 1
    end
    return SurfaceHoppingVariables(ComponentVector(v=v, r=r, σreal=σ, σimag=zero(σ)), state)
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

function evaluate_hopping_probability!(sim::RingPolymerSimulation{<:FSSH}, u, dt)
    v = RingPolymers.get_centroid(DynamicsUtils.get_velocities(u))
    σ = DynamicsUtils.get_quantum_subsystem(u)
    s = sim.method.state
    d = sim.calculator.centroid_nonadiabatic_coupling

    fewest_switches_probability!(sim.method.hopping_probability, v, σ, s, d, dt)
end

function rescale_velocity!(sim::RingPolymerSimulation{<:FSSH}, u)::Bool
    old_state = sim.method.state
    new_state = sim.method.new_state
    velocity = DynamicsUtils.get_velocities(u)
    centroid_velocity = RingPolymers.get_centroid(DynamicsUtils.get_velocities(u))

    c = calculate_potential_energy_change(sim.calculator, new_state, old_state)
    a, b = evaluate_a_and_b(sim, centroid_velocity, new_state, old_state)
    discriminant = b.^2 .- 2a.*c

    any(discriminant .< 0) && return false

    root = sqrt.(discriminant)
    plus = (b .+ root) ./ a
    minus = (b .- root) ./ a 
    velocity_rescale = sum(abs.(plus)) < sum(abs.(minus)) ? plus : minus
    perform_rescaling!(sim, velocity, velocity_rescale, new_state, old_state)

    return true
end

function evaluate_a_and_b(sim::RingPolymerSimulation{<:SurfaceHopping}, velocity, new_state, old_state)
    a = zeros(length(sim.atoms))
    b = zero(a)
    @views for i in range(sim.atoms)
        coupling = [sim.calculator.centroid_nonadiabatic_coupling[j,i][new_state, old_state] for j=1:ndofs(sim)]
        a[i] = dot(coupling, coupling) / sim.atoms.masses[i]
        b[i] = dot(velocity[:,i], coupling)
    end
    return (a, b)
end

function perform_rescaling!(sim::RingPolymerSimulation{<:SurfaceHopping}, velocity, velocity_rescale, new_state, old_state)
    for i in range(sim.atoms)
        coupling = [sim.calculator.centroid_nonadiabatic_coupling[j,i][new_state, old_state] for j=1:ndofs(sim)]
        velocity[:,i,:] .-= velocity_rescale[i] .* coupling ./ sim.atoms.masses[i]
    end
    return nothing
end

function calculate_potential_energy_change(calc::RingPolymerDiabaticCalculator, new_state::Integer, current_state::Integer)
    return calc.centroid_eigenvalues[new_state] - calc.centroid_eigenvalues[current_state]
end

function Estimators.diabatic_population(sim::RingPolymerSimulation{<:FSSH}, u)
    Calculators.evaluate_centroid_potential!(sim.calculator, DynamicsUtils.get_positions(u))
    U = eigvecs(sim.calculator.centroid_potential)

    σ = copy(DynamicsUtils.get_quantum_subsystem(u).re)
    σ[diagind(σ)] .= 0
    σ[u.state, u.state] = 1

    return diag(U * σ * U')
end

function DynamicsUtils.classical_hamiltonian(sim::RingPolymerSimulation{<:FSSH}, u)
    kinetic = DynamicsUtils.classical_kinetic_energy(sim, DynamicsUtils.get_velocities(u))
    spring = RingPolymers.get_spring_energy(sim.beads, sim.atoms.masses, DynamicsUtils.get_positions(u))

    Calculators.evaluate_potential!(sim.calculator, DynamicsUtils.get_positions(u))
    Calculators.eigen!(sim.calculator)
    potential = sum(bead[u.state] for bead in sim.calculator.eigenvalues)
    return kinetic + spring + potential
end
