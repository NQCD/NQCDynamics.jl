using StatsBase: mean
using ComponentArrays: ComponentVector
using NQCDistributions: ElectronicDistribution

export Ehrenfest

"""
    Ehrenfest{T} <: AbstractEhrenfest

Ehrenfest molecular dynamics. Classical molecular dynamics where the force is derived
by averaging contributions from multiple electronic states.

```jldoctest
Simulation{Ehrenfest}(Atoms(:H), DoubleWell())

# output

Simulation{Ehrenfest{Float64}}:
  Atoms{Float64}([:H], [1], [1837.4715941070515])
  DoubleWell{Int64, Int64, Int64, Int64}
  mass: Int64 1
  ω: Int64 1
  γ: Int64 1
  Δ: Int64 1
```
"""
struct Ehrenfest{T} <: AbstractEhrenfest
    density_propagator::Matrix{Complex{T}}
    tmp_complex_matrix::Matrix{Complex{T}}
    tmp_complex_matrix2::Matrix{Complex{T}}
    function Ehrenfest{T}(n_states::Integer) where {T}
        density_propagator = zeros(T, n_states, n_states)
        tmp_complex_matrix = zeros(Complex{T}, n_states, n_states)
        tmp_complex_matrix2 = zeros(Complex{T}, n_states, n_states)
        new{T}(density_propagator, tmp_complex_matrix, tmp_complex_matrix2)
    end
end

function Simulation{Ehrenfest}(atoms::Atoms{T}, model::Model; kwargs...) where {T}
    Simulation(atoms, model, Ehrenfest{T}(NQCModels.nstates(model)); kwargs...)
end

function DynamicsMethods.DynamicsVariables(
    sim::AbstractSimulation{<:AbstractEhrenfest}, v, r, electronic::ElectronicDistribution
)
    σ = DynamicsUtils.initialise_adiabatic_density_matrix(electronic, sim.cache, r)
    return ComponentVector(v=v, r=r, σreal=σ, σimag=zero(σ))
end

function DynamicsUtils.acceleration!(dv, v, r, sim::Simulation{<:Ehrenfest}, t, σ)
    fill!(dv, zero(eltype(dv)))
    NQCModels.state_independent_derivative!(sim.cache.model, dv, r)
    LinearAlgebra.lmul!(-1, dv)

    NQCCalculators.update_cache!(sim.cache, r) # Ensure adiabatic derivative is updated. 
    adiabatic_derivative = NQCCalculators.get_adiabatic_derivative(sim.cache, r)
    @inbounds for i in mobileatoms(sim)
        for j in dofs(sim)
            for m in eachstate(sim)
                for n in eachstate(sim)
                    dv[j,i] -= adiabatic_derivative[j,i][n,m] * real(σ[n,m])
                end
            end
        end
    end
    DynamicsUtils.divide_by_mass!(dv, sim.atoms.masses)
    return nothing
end

function Estimators.adiabatic_population(::AbstractSimulation{<:Ehrenfest}, u)
    σ = DynamicsUtils.get_quantum_subsystem(u)
    return real.(diag(σ))
end

function Estimators.diabatic_population(sim::AbstractSimulation{<:AbstractEhrenfest}, u)
    r = DynamicsUtils.get_positions(u)
    NQCDynamics.NQCCalculators.update_cache!(sim.cache, r) # Ensure eigenvecs needed for transformation are updated. 
    U = DynamicsUtils.evaluate_transformation(sim.cache, r)

    σ = DynamicsUtils.get_quantum_subsystem(u)

    return real.(diag(U * σ * U'))
end

function DynamicsUtils.classical_potential_energy(sim::Simulation{<:Ehrenfest}, u)
    NQCDynamics.NQCCalculators.update_cache!(sim.cache, DynamicsUtils.get_positions(u)) # Ensure eigen are updated
    eigs = NQCCalculators.get_eigen(sim.cache, DynamicsUtils.get_positions(u))

    potential = NQCModels.state_independent_potential(sim.cache.model, DynamicsUtils.get_positions(u))
    σ = DynamicsUtils.get_quantum_subsystem(u)
    for i in eachindex(eigs.values)
        potential += real(σ[i,i]) * eigs.values[i]
    end
    return potential
end
