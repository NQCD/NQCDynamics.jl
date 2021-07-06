using .Calculators: RingPolymerDiabaticCalculator

function acceleration!(dv, v, r, sim::RingPolymerSimulation{<:FSSH}, t, state)
    for i in axes(dv, 3)
        for j in axes(dv, 2)
            for k in axes(dv, 1)
                dv[k,j,i] = -sim.calculator.adiabatic_derivative[k,j,i][state, state] / sim.atoms.masses[j]
            end
        end
    end
    apply_interbead_coupling!(dv, r, sim)
    return nothing
end

function set_quantum_derivative!(dσ, v, σ, sim::RingPolymerSimulation{<:FSSH})
    V = sim.method.density_propagator

    V .= diagm(sum(sim.calculator.eigenvalues))
    for I in eachindex(v)
        @. V -= im * v[I] * sim.calculator.nonadiabatic_coupling[I]
    end
    V ./= length(sim.beads)
    mul!(sim.calculator.tmp_mat_complex1, V, σ)
    mul!(sim.calculator.tmp_mat_complex2, σ, V)
    @. dσ = -im * (sim.calculator.tmp_mat_complex1 - sim.calculator.tmp_mat_complex2)
    return nothing
end

function evaluate_hopping_probability!(sim::RingPolymerSimulation{<:FSSH}, u, dt)
    v = get_velocities(u)
    σ = get_quantum_subsystem(u)
    s = u.state
    d = sim.calculator.nonadiabatic_coupling

    sim.method.hopping_probability .= 0 # Set all entries to 0
    for m=1:sim.calculator.model.n_states
        if m != s
            for I in eachindex(v)
                sim.method.hopping_probability[m] += 2v[I]*real(σ[m,s]/σ[s,s])*d[I][s,m] * dt
            end
        end
    end
    sim.method.hopping_probability ./= length(sim.beads)
    clamp!(sim.method.hopping_probability, 0, 1)
    cumsum!(sim.method.hopping_probability, sim.method.hopping_probability)
    return nothing
end

function rescale_velocity!(sim::RingPolymerSimulation{<:FSSH}, u)::Bool
    old_state = u.state
    new_state = sim.method.new_state
    velocity = get_velocities(u)
    
    c = calculate_potential_energy_change(sim.calculator, new_state, old_state)
    a, b = evaluate_a_and_b(sim, velocity, new_state, old_state)
    discriminant = zero(a)
    for i in axes(discriminant, 2)
        for j in axes(discriminant, 1)
            discriminant[j,i] = b[j,i]^2 - 2a[j,i] * c[i]
        end
    end

    any(discriminant .< 0) && return false

    root = sqrt.(discriminant)
    plus = (b .+ root) ./ a
    minus = (b .+ root) ./ a
    velocity_rescale = sum(abs.(plus)) < sum(abs.(minus)) ? plus : minus

    perform_rescaling!(sim, velocity, velocity_rescale, new_state, old_state)

    return true
end

function evaluate_a_and_b(sim::RingPolymerSimulation{<:FSSH}, velocity, new_state, old_state)
    a = zeros(length(sim.atoms), length(sim.beads))
    b = zero(a)
    @views for i in range(sim.beads)
        for j in range(sim.atoms)
            coupling = [sim.calculator.nonadiabatic_coupling[k,j,i][new_state, old_state] for k=1:sim.DoFs]
            a[j,i] = dot(coupling, coupling) / sim.atoms.masses[j]
            b[j,i] = dot(velocity[:,j,i], coupling)
        end
    end
    return (a, b)
end

function perform_rescaling!(sim::RingPolymerSimulation{<:FSSH}, velocity, velocity_rescale, new_state, old_state)
    for i=1:length(sim.beads)
        for j=1:length(sim.atoms)
            coupling = [sim.calculator.nonadiabatic_coupling[k,j,i][new_state, old_state] for k=1:sim.DoFs]
            velocity[:,j,i] .-= velocity_rescale[j,i] .* coupling ./ sim.atoms.masses[j]
        end
    end
    return nothing
end

function calculate_potential_energy_change(calc::RingPolymerDiabaticCalculator, new_state::Integer, current_state::Integer)
    return [eigs[new_state] - eigs[current_state] for eigs in calc.eigenvalues]
end

function get_diabatic_population(sim::RingPolymerSimulation{<:FSSH}, u)
    Calculators.evaluate_centroid_potential!(sim.calculator, get_positions(u))
    U = eigvecs(sim.calculator.potential[1])

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
