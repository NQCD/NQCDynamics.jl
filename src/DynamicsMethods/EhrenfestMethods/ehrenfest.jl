using StatsBase: mean
using .Calculators: DiabaticCalculator, RingPolymerDiabaticCalculator
using ComponentArrays: ComponentVector

export Ehrenfest

struct Ehrenfest{T} <: AbstractEhrenfest
    density_propagator::Matrix{Complex{T}}
    function Ehrenfest{T}(n_states::Integer) where {T}
        density_propagator = zeros(n_states, n_states)
        new{T}(density_propagator)
    end
end

function Simulation{Ehrenfest}(atoms::Atoms{S,T}, model::Model; kwargs...) where {S,T}
    Simulation(atoms, model, Ehrenfest{T}(NonadiabaticModels.nstates(model)); kwargs...)
end

function DynamicsMethods.DynamicsVariables(sim::Simulation{<:AbstractEhrenfest}, v, r, state::Integer; type=:diabatic)
    n_states = NonadiabaticModels.nstates(sim.calculator.model)
    if type == :diabatic
        Calculators.evaluate_potential!(sim.calculator, r)
        Calculators.eigen!(sim.calculator)
        U = sim.calculator.eigenvectors

        diabatic_density = zeros(n_states, n_states)
        diabatic_density[state, state] = 1
        σ = U' * diabatic_density * U
    else
        σ = zeros(n_states, n_states)
        σ[state, state] = 1
    end
    return ComponentVector(v=v, r=r, σreal=σ, σimag=zero(σ))
end

function acceleration!(dv, v, r, sim::AbstractSimulation{<:Ehrenfest}, t, σ)
    dv .= zero(eltype(dv))
    for I in eachindex(dv)
        for J in eachindex(σ)
            dv[I] -= sim.calculator.adiabatic_derivative[I][J] * real(σ[J])
        end
    end
    DynamicsUtils.divide_by_mass!(dv, sim.atoms.masses)
    return nothing
end

function Estimators.adiabatic_population(::AbstractSimulation{<:Ehrenfest}, u)
    σ = DynamicsUtils.get_quantum_subsystem(u)
    return real.(diag(σ))
end

function Estimators.diabatic_population(sim::Simulation{<:Ehrenfest}, u)
    Calculators.evaluate_potential!(sim.calculator, DynamicsUtils.get_positions(u))
    U = eigvecs(sim.calculator.potential)

    σ = DynamicsUtils.get_quantum_subsystem(u)

    return real.(diag(U * σ * U'))
end

function NonadiabaticMolecularDynamics.evaluate_hamiltonian(sim::Simulation{<:Ehrenfest}, u)
    k = NonadiabaticMolecularDynamics.evaluate_kinetic_energy(sim.atoms.masses, DynamicsUtils.get_velocities(u))
    Calculators.evaluate_potential!(sim.calculator, DynamicsUtils.get_positions(u))
    Calculators.eigen!(sim.calculator)
    p = sum(diag(DynamicsUtils.get_quantum_subsystem(u)) .* sim.calculator.eigenvalues)
    return k + p
end
