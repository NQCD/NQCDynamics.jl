using LinearAlgebra: lmul!

"""
    CMM2{T} <: DynamicsMethods.Method

Classical mapping model 2 from Xin He and Jian Liu (2019).
"""
struct CMM2{T} <: DynamicsMethods.Method
    temp_q::Vector{T}
    temp_p::Vector{T}
    function CMM2{T}(n_states::Integer) where {T}
        new{T}(zeros(n_states), zeros(n_states))
    end
end

function NonadiabaticMolecularDynamics.Simulation{CMM2}(atoms::Atoms{S,T}, model::Model; kwargs...) where {S,T}
    NonadiabaticMolecularDynamics.Simulation(atoms, model, CMM2{T}(NonadiabaticModels.nstates(model)); kwargs...)
end

function DynamicsMethods.DynamicsVariables(sim::Simulation{<:CMM2}, v, r, state::Integer; type=:diabatic)

    qmap = zeros(NonadiabaticModels.nstates(sim))
    pmap = zeros(NonadiabaticModels.nstates(sim))

    θ = rand() * 2π
    qmap[state] = cos(θ) * sqrt(2)
    pmap[state] = sin(θ) * sqrt(2)

    ComponentVector(v=v, r=r, pmap=pmap, qmap=qmap)
end

function DynamicsMethods.motion!(du, u, sim::Simulation{<:CMM2}, t)
    dr = DynamicsUtils.get_positions(du)
    dv = DynamicsUtils.get_velocities(du)
    r = DynamicsUtils.get_positions(u)
    v = DynamicsUtils.get_velocities(u)

    DynamicsUtils.velocity!(dr, v, r, sim, t)
    acceleration!(dv, u, sim)
    set_mapping_force!(du, u, sim)
end

function acceleration!(dv, u, sim::Simulation{<:CMM2})

    Calculators.evaluate_derivative!(sim.calculator, DynamicsUtils.get_positions(u))
    qmap = get_mapping_positions(u)
	pmap = get_mapping_momenta(u)
    for I in eachindex(dv)
        D = sim.calculator.derivative[I]
        mul!(sim.method.temp_q, D, qmap)
        mul!(sim.method.temp_p, D, pmap)
        dv[I] = dot(qmap, sim.method.temp_q)
        dv[I] += dot(pmap, sim.method.temp_p)
    end
    lmul!(-1/2, dv)
    DynamicsUtils.divide_by_mass!(dv, sim.atoms.masses)
end

function set_mapping_force!(du, u, sim::Simulation{<:CMM2})
    Calculators.evaluate_potential!(sim.calculator, DynamicsUtils.get_positions(u))
    V = sim.calculator.potential
    mul!(get_mapping_positions(du), V, get_mapping_momenta(u))
    mul!(get_mapping_momenta(du), V, get_mapping_positions(u))
    lmul!(-1, get_mapping_momenta(du))
end

function Estimators.diabatic_population(sim::Simulation{<:CMM2}, u)
    qmap = get_mapping_positions(u)
    pmap = get_mapping_momenta(u)
    out = zero(qmap)
    for i=1:NonadiabaticModels.nstates(sim)
        out[i] = (qmap[i]^2 + pmap[i]^2) / 2
    end
    out
end

function DynamicsUtils.classical_hamiltonian(sim::Simulation{<:CMM2}, u)
    r = DynamicsUtils.get_positions(u)
    v = DynamicsUtils.get_velocities(u)

    kinetic = DynamicsUtils.classical_kinetic_energy(sim, v)

    Calculators.evaluate_potential!(sim.calculator, r)
    V = sim.calculator.potential
    qmap = get_mapping_positions(u)
    pmap = get_mapping_momenta(u)
    potential = (pmap'V*pmap + qmap'V*qmap) / 2

    H = kinetic + potential
    return H
end
