using LinearAlgebra: lmul!, norm
using Parameters: Parameters
using NQCDynamics.NonadiabaticDistributions: ElectronicDistribution

"""
    eCMM{T} <: DynamicsMethods.Method

# References

- [HeGong2021](@cite)
- [HeWu2021](@cite)
"""
struct eCMM{T} <: DynamicsMethods.Method
    γ::T
    temp_q::Vector{T}
    temp_p::Vector{T}
    function eCMM{T}(n_states::Integer, γ) where {T}
        new{T}(γ, zeros(n_states), zeros(n_states))
    end
end

function NQCDynamics.Simulation{eCMM}(atoms::Atoms{T}, model::Model; γ=0, kwargs...) where {T}
    NQCDynamics.Simulation(atoms, model, eCMM{T}(NQCModels.nstates(model), γ); kwargs...)
end

function DynamicsMethods.DynamicsVariables(sim::Simulation{<:eCMM}, v, r, ::ElectronicDistribution)

    F = NQCModels.nstates(sim)
    radius = sqrt(2 + 2F*sim.method.γ)
    points = generate_random_points_on_nsphere(2F, radius)
    qmap = points[begin:div(F,2)+1]
    pmap = points[div(F,2)+2:end]

    ComponentVector(v=v, r=r, pmap=pmap, qmap=qmap)
end

function generate_random_points_on_nsphere(n, radius)
    randomnormal = randn(n)
    return randomnormal ./ norm(randomnormal) .* radius
end

function DynamicsMethods.motion!(du, u, sim::AbstractSimulation{<:eCMM}, t)
    dr = DynamicsUtils.get_positions(du)
    dv = DynamicsUtils.get_velocities(du)
    r = DynamicsUtils.get_positions(u)
    v = DynamicsUtils.get_velocities(u)

    DynamicsUtils.velocity!(dr, v, r, sim, t)
    acceleration!(dv, u, sim)
    set_mapping_force!(du, u, sim)
end

function acceleration!(dv, u, sim::Simulation{<:eCMM})

    Parameters.@unpack γ, temp_q, temp_p = sim.method

    Calculators.evaluate_derivative!(sim.calculator, DynamicsUtils.get_positions(u))
    qmap = get_mapping_positions(u)
	pmap = get_mapping_momenta(u)
    for I in eachindex(dv)
        D = sim.calculator.derivative[I]
        mul!(temp_q, D, qmap)
        mul!(temp_p, D, pmap)
        dv[I] = dot(qmap, temp_q)
        dv[I] += dot(pmap, temp_p)
        dv[I] -= tr(D) * 2γ
    end
    lmul!(-1/2, dv)
    DynamicsUtils.divide_by_mass!(dv, sim.atoms.masses)
end

function set_mapping_force!(du, u, sim::Simulation{<:eCMM})
    Calculators.evaluate_potential!(sim.calculator, DynamicsUtils.get_positions(u))
    V = sim.calculator.potential
    mul!(get_mapping_positions(du), V, get_mapping_momenta(u))
    mul!(get_mapping_momenta(du), V, get_mapping_positions(u))
    lmul!(-1, get_mapping_momenta(du))
end

function DynamicsUtils.classical_hamiltonian(sim::Simulation{<:eCMM}, u)
    r = DynamicsUtils.get_positions(u)
    v = DynamicsUtils.get_velocities(u)

    kinetic = DynamicsUtils.classical_kinetic_energy(sim, v)

    Calculators.evaluate_potential!(sim.calculator, r)
    V = sim.calculator.potential
    qmap = get_mapping_positions(u)
    pmap = get_mapping_momenta(u)
    potential = (pmap'V*pmap + qmap'V*qmap - 2sim.method.γ*tr(V)) / 2

    H = kinetic + potential
    return H
end

mapping_kernel(qmap, pmap, γ) = (qmap.^2 .+ pmap.^2)./2 .- γ

function inverse_mapping_kernel(qmap, pmap, γ)
    F = length(qmap)
    prefactor = (1 + F) / (1 + F*γ)^2
    subtractor = (1 - γ) / (1 + F*γ)
    return prefactor .* (qmap.^2 .+ pmap.^2)./2 .- subtractor
end

function Estimators.diabatic_population(sim::Simulation{<:eCMM}, u)
    qmap = get_mapping_positions(u)
    pmap = get_mapping_momenta(u)
    inverse_mapping_kernel(qmap, pmap, sim.method.γ)
end

function Estimators.initial_diabatic_population(sim::Simulation{<:eCMM}, u)
    qmap = get_mapping_positions(u)
    pmap = get_mapping_momenta(u)
    mapping_kernel(qmap, pmap, sim.method.γ)
end

function TimeCorrelationFunctions.evaluate_normalisation(
    sim::AbstractSimulation{<:eCMM},
    ::TimeCorrelationFunctions.PopulationCorrelationFunction
)
    NQCModels.nstates(sim)
end

