using LinearAlgebra: eigvecs, diag, diagind, dot
using StatsBase: sample, Weights
using RingPolymerArrays: get_centroid

using NQCDynamics.Calculators: RingPolymerDiabaticCalculator
using NQCDynamics: RingPolymerSimulation, RingPolymers

function RingPolymerSimulation{FSSH}(atoms::Atoms{T}, model::Model, n_beads::Integer; rescaling=:standard, kwargs...) where {T}
    RingPolymerSimulation(atoms, model, FSSH{T}(NQCModels.nstates(model), rescaling), n_beads; kwargs...)
end

function DynamicsMethods.motion!(du, u, sim::RingPolymerSimulation{<:SurfaceHopping}, t)
    dr = DynamicsUtils.get_positions(du)
    dv = DynamicsUtils.get_velocities(du)
    dσ = DynamicsUtils.get_quantum_subsystem(du)

    r = DynamicsUtils.get_positions(u)
    v = DynamicsUtils.get_velocities(u)
    σ = DynamicsUtils.get_quantum_subsystem(u)

    set_state!(u, sim.method.state) # Make sure the state variables match, 

    DynamicsUtils.velocity!(dr, v, r, sim, t)
    Calculators.update_electronics!(sim.calculator, r)
    DynamicsUtils.acceleration!(dv, v, r, sim, t, sim.method.state)
    DynamicsUtils.apply_interbead_coupling!(dv, r, sim)
    DynamicsUtils.set_quantum_derivative!(dσ, u, sim)
end

function DynamicsUtils.get_hopping_nonadiabatic_coupling(sim::RingPolymerSimulation, r::AbstractArray{T,3}) where {T}
    return Calculators.get_centroid_nonadiabatic_coupling(sim.calculator, r)
end

function DynamicsUtils.get_hopping_velocity(::RingPolymerSimulation, v::AbstractArray{T,3}) where {T}
    return get_centroid(v)
end

function DynamicsUtils.get_hopping_eigenvalues(sim::RingPolymerSimulation, r::AbstractArray{T,3}) where {T}
    return Calculators.get_centroid_eigen(sim.calculator, r).values
end

function perform_rescaling!(
    sim::RingPolymerSimulation{<:SurfaceHopping}, velocity, γ, d
)
    for I in CartesianIndices(d)
        velocity[I,:] .-= γ * d[I] / masses(sim, I)
    end
    return nothing
end

function DynamicsUtils.classical_potential_energy(sim::RingPolymerSimulation{<:FSSH}, u)
    all_eigs = Calculators.get_eigen(sim.calculator, DynamicsUtils.get_positions(u))
    potential = sum(eigs.values[u.state] for eigs in all_eigs)
    return potential
end

function DynamicsUtils.centroid_classical_potential_energy(sim::RingPolymerSimulation{<:FSSH}, u)
    centroid_eigs = Calculators.get_centroid_eigen(sim.calculator, DynamicsUtils.get_positions(u))
    potential = centroid_eigs.values[u.state]
    return potential
end
