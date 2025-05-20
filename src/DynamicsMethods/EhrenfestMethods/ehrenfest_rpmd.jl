using NQCDynamics: RingPolymers
using NQCCalculators
function RingPolymerSimulation{Ehrenfest}(atoms::Atoms{T}, model::Model, n_beads::Integer; kwargs...) where {T}
    RingPolymerSimulation(atoms, model, Ehrenfest{T}(NQCModels.nstates(model)), n_beads; kwargs...)
end

function DynamicsMethods.motion!(du, u, sim::RingPolymerSimulation{<:AbstractEhrenfest}, t)
    dr = DynamicsUtils.get_positions(du)
    dv = DynamicsUtils.get_velocities(du)
    dσ = DynamicsUtils.get_quantum_subsystem(du)

    r = DynamicsUtils.get_positions(u)
    v = DynamicsUtils.get_velocities(u)
    σ = DynamicsUtils.get_quantum_subsystem(u)

    DynamicsUtils.velocity!(dr, v, r, sim, t)
    NQCCalculators.update_electronics!(sim.cache, r)
    DynamicsUtils.acceleration!(dv, v, r, sim, t, σ)
    DynamicsUtils.apply_interbead_coupling!(dv, r, sim)
    DynamicsUtils.set_quantum_derivative!(dσ, u, sim)
end

function DynamicsUtils.acceleration!(dv, v, r, sim::RingPolymerSimulation{<:Ehrenfest}, t, σ)
    fill!(dv, zero(eltype(dv)))
    NQCModels.state_independent_derivative!(sim.cache.model, dv, r)
    LinearAlgebra.lmul!(-1, dv)

    adiabatic_derivative = NQCCalculators.get_adiabatic_derivative(sim.cache, r)
    for b in axes(dv, 3)
        for i in mobileatoms(sim)
            for j in dofs(sim)
                for m in eachstate(sim)
                    for n in eachstate(sim)
                        dv[j,i,b] -= adiabatic_derivative[j,i,b][n,m] * real(σ[n,m])
                    end
                end
            end
        end
    end
    DynamicsUtils.divide_by_mass!(dv, sim.atoms.masses)

    return nothing
end

function DynamicsUtils.classical_potential_energy(sim::RingPolymerSimulation{<:Ehrenfest}, u)
    r = DynamicsUtils.get_positions(u)
    all_eigs = NQCCalculators.get_eigen(sim.cache, r)
    population = Estimators.adiabatic_population(sim, u)
    potential = sum(dot(population, eigs.values) for eigs in all_eigs)
    return potential
end
