using LinearAlgebra: lmul!
using Parameters: Parameters

"""
    eCMM{T} <: DynamicsMethods.Method

# References

[He2021a](@cite)
[He2021b](@cite)
"""
struct eCMM{T} <: DynamicsMethods.Method
    γ::T
    temp_q::Vector{T}
    temp_p::Vector{T}
    function eCMM{T}(n_states::Integer, γ) where {T}
        new{T}(γ, zeros(n_states), zeros(n_states))
    end
end

function NonadiabaticMolecularDynamics.Simulation{eCMM}(atoms::Atoms{S,T}, model::Model; γ=0, kwargs...) where {S,T}
    NonadiabaticMolecularDynamics.Simulation(atoms, model, eCMM{T}(NonadiabaticModels.nstates(model), γ); kwargs...)
end

function DynamicsMethods.DynamicsVariables(sim::Simulation{<:eCMM}, v, r, state::Integer; type=:diabatic)

    qmap = zeros(NonadiabaticModels.nstates(sim))
    pmap = zeros(NonadiabaticModels.nstates(sim))

    θ = rand() * 2π
    radius = sqrt(2 + 2NonadiabaticModels.nstates(sim)*sim.method.γ)
    qmap[state] = cos(θ) * radius
    pmap[state] = sin(θ) * radius

    ComponentVector(v=v, r=r, pmap=pmap, qmap=qmap)
end

function DynamicsMethods.motion!(du, u, sim::Simulation{<:eCMM}, t)
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

function Estimators.diabatic_population(sim::Simulation{<:eCMM}, u)
    qmap = get_mapping_positions(u)
    pmap = get_mapping_momenta(u)
    F = NonadiabaticModels.nstates(sim)
    γ = sim.method.γ

    prefactor = (1 + F)/(2*(1 + F*γ)^2)
    subtractor = (1 - γ) / (1 + F*γ)
    out = zero(qmap)
    for i=1:F
        out[i] = prefactor * (qmap[i]^2 + pmap[i]^2) - subtractor
    end
    out
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
