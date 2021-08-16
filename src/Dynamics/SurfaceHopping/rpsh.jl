using .Calculators: RingPolymerDiabaticCalculator

function DynamicsVariables(sim::RingPolymerSimulation{<:SurfaceHopping}, v, r, state::Integer; type=:diabatic)
    n_states = sim.calculator.model.n_states
    if type == :diabatic
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

function acceleration!(dv, v, r, sim::RingPolymerSimulation{<:FSSH}, t, state)
    for I in eachindex(dv)
        dv[I] = -sim.calculator.adiabatic_derivative[I][state, state]
    end
    divide_by_mass!(dv, sim.atoms.masses)
    apply_interbead_coupling!(dv, r, sim)
    return nothing
end

function evaluate_hopping_probability!(sim::RingPolymerSimulation{<:FSSH}, u, dt)
    v = get_centroid(get_velocities(u))
    σ = get_quantum_subsystem(u)
    s = sim.method.state
    d = sim.calculator.centroid_nonadiabatic_coupling

    fewest_switches_probability!(sim.method.hopping_probability, v, σ, s, d, dt)
end

function rescale_velocity!(sim::RingPolymerSimulation{<:FSSH}, u)::Bool
    old_state = sim.method.state
    new_state = sim.method.new_state
    velocity = get_velocities(u)
    centroid_velocity = get_centroid(get_velocities(u))

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
        coupling = [sim.calculator.centroid_nonadiabatic_coupling[j,i][new_state, old_state] for j=1:sim.DoFs]
        a[i] = dot(coupling, coupling) / sim.atoms.masses[i]
        b[i] = dot(velocity[:,i], coupling)
    end
    return (a, b)
end

function perform_rescaling!(sim::RingPolymerSimulation{<:SurfaceHopping}, velocity, velocity_rescale, new_state, old_state)
    for i in range(sim.atoms)
        coupling = [sim.calculator.centroid_nonadiabatic_coupling[j,i][new_state, old_state] for j=1:sim.DoFs]
        velocity[:,i,:] .-= velocity_rescale[i] .* coupling ./ sim.atoms.masses[i]
    end
    return nothing
end

function calculate_potential_energy_change(calc::RingPolymerDiabaticCalculator, new_state::Integer, current_state::Integer)
    return calc.centroid_eigenvalues[new_state] - calc.centroid_eigenvalues[current_state]
end

function get_diabatic_population(sim::RingPolymerSimulation{<:FSSH}, u)
    Calculators.evaluate_centroid_potential!(sim.calculator, get_positions(u))
    U = eigvecs(sim.calculator.centroid_potential)

    σ = copy(get_quantum_subsystem(u))
    σ[diagind(σ)] .= 0
    σ[u.state, u.state] = 1

    return real.(diag(U * σ * U'))
end

function NonadiabaticMolecularDynamics.evaluate_hamiltonian(sim::RingPolymerSimulation{<:FSSH}, u)
    k = evaluate_kinetic_energy(sim.atoms.masses, get_velocities(u))
    Calculators.evaluate_potential!(sim.calculator, get_positions(u))
    Calculators.eigen!(sim.calculator)
    p = sum([bead[u.state] for bead in sim.calculator.eigenvalues])
    return k + p + get_spring_energy(sim, get_positions(u))
end
