
function acceleration!(dv, v, r, sim::RingPolymerSimulation{<:FSSH}, t; state=1)
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

function set_density_matrix_derivative!(dσ, v, σ, sim::RingPolymerSimulation{<:FSSH})
    V = sim.method.density_propagator

    V .= diagm(mean(sim.calculator.eigenvalues))
    ivd = mean(im .* v .* sim.calculator.nonadiabatic_coupling; dims=3)
    for I in eachindex(ivd)
        V .-= ivd[I]
    end
    dσ .= -im*(V*σ - σ*V)
    return nothing
end

function evaluate_hopping_probability!(sim::RingPolymerSimulation{<:FSSH}, u, dt)
    v = get_velocities(u)
    σ = get_density_matrix(u)
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
    return nothing
end

function evaluate_a_and_b(sim::RingPolymerSimulation{<:FSSH}, velocity::RingPolymerArray, new_state, old_state)
    a = zeros(length(sim.atoms), length(sim.beads))
    b = zero(a)
    @views for i in range(sim.beads)
        for j in range(sim.atoms)
            coupling = [sim.calculator.nonadiabatic_coupling[k,j,i][new_state, old_state] for k=1:sim.DoFs]
            a[j,i] = coupling'coupling / sim.atoms.masses[j]
            b[j,i] = velocity[:,j,i]'coupling
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
    return mean([eigs[new_state] - eigs[current_state] for eigs in calc.eigenvalues])
end

function get_population(sim::RingPolymerSimulation{<:FSSH}, u)
    Models.potential!(sim.calculator.model, sim.calculator.potential[1], dropdims(mean(get_positions(u); dims=3), dims=3))
    vals, U = eigen!(sim.calculator.potential[1])

    σ = copy(get_density_matrix(u))
    σ[diagind(σ)] .= 0
    σ[u.state, u.state] = 1

    return real.(diag(U * σ * U'))
end