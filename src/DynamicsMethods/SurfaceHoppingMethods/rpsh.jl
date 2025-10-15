using LinearAlgebra: eigvecs, diag, diagind, dot
using StatsBase: sample, Weights
using RingPolymerArrays: get_centroid

using NQCDynamics: RingPolymerSimulation, RingPolymers
using NQCCalculators

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

    set_state!(u, sim.method.state, sim) # Make sure the state variables match, 

    DynamicsUtils.velocity!(dr, v, r, sim, t)
    NQCCalculators.update_cache!(sim.cache, r)
    DynamicsUtils.acceleration!(dv, v, r, sim, t, sim.method.state)
    DynamicsUtils.apply_interbead_coupling!(dv, r, sim)
    DynamicsUtils.set_quantum_derivative!(dσ, u, sim)
end

function perform_rescaling!(
    sim::RingPolymerSimulation{<:SurfaceHopping}, velocity, γ, d
)
    for I in CartesianIndices(d)
        velocity[I,:] .-= γ * d[I] / masses(sim, I)
    end
    return nothing
end

function frustrated_hop_invert_velocity!(
    sim::RingPolymerSimulation{<:SurfaceHopping}, velocity, d
)
    #invert center of mass velocity of ring polymer
    dn = LinearAlgebra.normalize(d)
    vcom = sum(velocity, dims=3) / size(velocity,3)
    γ = dot(vcom[:,:,1],dn)
    for I in CartesianIndices(dn)
        velocity[I,:] .-= 2γ * dn[I]
    end
    return nothing
end

function DynamicsUtils.classical_potential_energy(sim::RingPolymerSimulation{<:FSSH}, u)
    NQCCalculators.update_cache!(sim.cache, DynamicsUtils.get_positions(u)) # Ensure eigen is populated
    all_eigs = NQCCalculators.get_eigen(sim.cache, DynamicsUtils.get_positions(u))
    potential = sum(eigs.values[convert(Int, u.state |> first)] for eigs in all_eigs)
    return potential
end

function DynamicsUtils.centroid_classical_potential_energy(sim::RingPolymerSimulation{<:FSSH}, u)
    NQCCalculators.update_cache!(sim.cache, DynamicsUtils.get_positions(u)) # Ensure eigen is populated
    centroid_eigs = NQCCalculators.get_centroid_eigen(sim.cache, DynamicsUtils.get_positions(u))
    potential = centroid_eigs.values[convert(Int, u.state |> first)]
    return potential
end
