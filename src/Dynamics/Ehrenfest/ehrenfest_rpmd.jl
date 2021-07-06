function motion!(du, u, sim::RingPolymerSimulation{<:AbstractEhrenfest}, t)
    dr = get_positions(du)
    dv = get_velocities(du)
    dσ = get_quantum_subsystem(du)
    r = get_positions(u)
    v = get_velocities(u)
    σ = get_quantum_subsystem(u)
    velocity!(dr, v, r, sim, t)
    Calculators.update_electronics!(sim.calculator, r)
    Calculators.update_centroid_electronics!(sim.calculator, r)
    acceleration!(dv, v, r, sim, t, σ)
    set_quantum_derivative!(dσ, get_centroid(v), σ, sim)
end

function acceleration!(dv, v, r, sim::RingPolymerSimulation{<:Ehrenfest}, t, σ)
    dv .= zero(eltype(dv))
    for I in eachindex(dv)
        for J in eachindex(σ)
            dv[I] -= sim.calculator.adiabatic_derivative[I][J] * real(σ[J])
        end
    end
    divide_by_mass!(dv, sim.atoms.masses)
    apply_interbead_coupling!(dv, r, sim)
    return nothing
end

function set_quantum_derivative_old!(dσ, v, σ, sim::RingPolymerSimulation{<:Ehrenfest})
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

function set_quantum_derivative!(dσ, v, σ, sim::RingPolymerSimulation{<:Ehrenfest})
    V = sim.method.density_propagator

    V .= diagm(sim.calculator.centroid_eigenvalues)
    for I in eachindex(v)
        @. V -= im * v[I] * sim.calculator.centroid_nonadiabatic_coupling[I]
    end
    #V ./= length(sim.beads)
    mul!(sim.calculator.tmp_mat_complex1, V, σ)
    mul!(sim.calculator.tmp_mat_complex2, σ, V)
    @. dσ = -im * (sim.calculator.tmp_mat_complex1 - sim.calculator.tmp_mat_complex2)
    return nothing
end

function get_diabatic_population(sim::RingPolymerSimulation{<:Ehrenfest}, u)
    Calculators.evaluate_centroid_potential!(sim.calculator, get_positions(u))
    #U = eigvecs(sim.calculator.potential[1])
    U = eigvecs(sim.calculator.centroid_potential)
    σ = get_quantum_subsystem(u)
    return real.(diag(U * σ * U'))
end

# function get_diabatic_population(sim::RingPolymerSimulation{<:Ehrenfest}, u)
#     Calculators.evaluate_potential!(sim.calculator, get_positions(u))
#     Calculators.eigen!(sim.calculator)
#     population = zeros(sim.calculator.model.n_states)
#     σ = get_quantum_subsystem(u)
#     for i=1:length(sim.beads)
#         U = sim.calculator.eigenvectors[i]
#         population .+= real.(diag(U * σ * U'))
#     end
#     population ./= length(sim.beads)
#     #U = eigvecs(sim.calculator.potential[1])

#     return population #real.(diag(U * σ * U'))
# end

function NonadiabaticMolecularDynamics.evaluate_hamiltonian(sim::RingPolymerSimulation{<:Ehrenfest}, u)
    k = evaluate_kinetic_energy(sim.atoms.masses, get_velocities(u))
    Calculators.evaluate_potential!(sim.calculator, get_positions(u))
    Calculators.eigen!(sim.calculator)

    population = get_adiabatic_population(sim, u)
    p = sum([dot(population, eigs) for eigs in sim.calculator.eigenvalues])
    return k + p
end
