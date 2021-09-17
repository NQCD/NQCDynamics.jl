
function RingPolymerSimulation{Ehrenfest}(atoms::Atoms{S,T}, model::Model, n_beads::Integer; kwargs...) where {S,T}
    RingPolymerSimulation(atoms, model, Ehrenfest{T}(NonadiabaticModels.nstates(model)), n_beads; kwargs...)
end

function DynamicsMethods.DynamicsVariables(sim::RingPolymerSimulation{<:AbstractEhrenfest}, v, r, state::Integer; type=:diabatic)
    n_states = NonadiabaticModels.nstates(sim.calculator.model)
    if type == :diabatic
        Calculators.evaluate_centroid_potential!(sim.calculator, r)
        U = eigvecs(sim.calculator.centroid_potential)

        diabatic_density = zeros(n_states, n_states)
        diabatic_density[state, state] = 1
        σ = U' * diabatic_density * U
    else
        σ = zeros(n_states, n_states)
        σ[state, state] = 1
    end
    return ComponentVector(v=v, r=r, σreal=σ, σimag=zero(σ))
end

function acceleration!(dv, v, r, sim::RingPolymerSimulation{<:Ehrenfest}, t, σ)
    dv .= zero(eltype(dv))
    for I in eachindex(dv)
        for J in eachindex(σ)
            dv[I] -= sim.calculator.adiabatic_derivative[I][J] * real(σ[J])
        end
    end
    DynamicsUtils.divide_by_mass!(dv, sim.atoms.masses)
    DynamicsUtils.apply_interbead_coupling!(dv, r, sim)
    return nothing
end

function Estimators.diabatic_population(sim::RingPolymerSimulation{<:Ehrenfest}, u)
    Calculators.evaluate_centroid_potential!(sim.calculator, DynamicsUtils.get_positions(u))
    U = eigvecs(sim.calculator.centroid_potential)
    σ = DynamicsUtils.get_quantum_subsystem(u)
    return real.(diag(U * σ * U'))
end

function NonadiabaticMolecularDynamics.evaluate_hamiltonian(sim::RingPolymerSimulation{<:Ehrenfest}, u)
    k = NonadiabaticMolecularDynamics.evaluate_kinetic_energy(sim.atoms.masses, DynamicsUtils.get_velocities(u))
    Calculators.evaluate_potential!(sim.calculator, DynamicsUtils.get_positions(u))
    Calculators.eigen!(sim.calculator)

    population = Estimators.adiabatic_population(sim, u)
    p = sum([dot(population, eigs) for eigs in sim.calculator.eigenvalues])
    return k + p
end
