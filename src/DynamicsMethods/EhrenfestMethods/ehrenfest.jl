using StatsBase: mean
using .NonadiabaticDistributions: NonadiabaticDistribution
using ComponentArrays: ComponentVector

export Ehrenfest

"""
    Ehrenfest{T} <: AbstractEhrenfest

Ehrenfest molecular dynamics. Classical molecular dynamics where the force is derived
by averaging contributions from multiple electronic states.

```jldoctest
Simulation{Ehrenfest}(Atoms(:H), DoubleWell())

# output

Simulation{Ehrenfest{Float64}}:
  Atoms{1, Float64}([:H], UInt8[0x01], [1837.4715941070515])
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
    function Ehrenfest{T}(n_states::Integer) where {T}
        density_propagator = zeros(T, n_states, n_states)
        tmp_complex_matrix = zeros(Complex{T}, n_states, n_states)
        new{T}(density_propagator, tmp_complex_matrix)
    end
end

function Simulation{Ehrenfest}(atoms::Atoms{S,T}, model::Model; kwargs...) where {S,T}
    Simulation(atoms, model, Ehrenfest{T}(NQCModels.nstates(model)); kwargs...)
end

function DynamicsMethods.DynamicsVariables(
    sim::AbstractSimulation{<:AbstractEhrenfest}, v, r, electronic::NonadiabaticDistribution
)
    σ = NonadiabaticDistributions.initialise_adiabatic_density_matrix(electronic, sim.calculator, r)
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

function Estimators.diabatic_population(sim::AbstractSimulation{<:AbstractEhrenfest}, u)
    r = DynamicsUtils.get_positions(u)
    U = NonadiabaticDistributions.evaluate_transformation(sim.calculator, r)

    σ = DynamicsUtils.get_quantum_subsystem(u)

    return real.(diag(U * σ * U'))
end

function DynamicsUtils.classical_hamiltonian(sim::Simulation{<:Ehrenfest}, u)
    kinetic = DynamicsUtils.classical_kinetic_energy(sim, DynamicsUtils.get_velocities(u))

    eigs = Calculators.get_eigen(sim.calculator, DynamicsUtils.get_positions(u))
    potential = sum(diag(DynamicsUtils.get_quantum_subsystem(u)) .* eigs.values)

    return kinetic + potential
end
