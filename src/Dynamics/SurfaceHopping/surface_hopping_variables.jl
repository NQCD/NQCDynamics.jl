using StatsBase: sample, Weights

mutable struct SurfaceHoppingVariables{T,A,Axes,S} <: DEDataVector{T}
    x::ComponentVector{T,A,Axes}
    state::S
end

function DynamicsVariables(sim::Simulation{<:SurfaceHopping}, v, r, state::Integer; type=:diabatic)
    n_states = sim.calculator.model.n_states
    if type == :diabatic
        Calculators.evaluate_potential!(sim.calculator, r)
        Calculators.eigen!(sim.calculator)
        U = sim.calculator.eigenvectors

        diabatic_density = zeros(n_states, n_states)
        diabatic_density[state, state] = 1
        σ = U' * diabatic_density * U
        state = sample(Weights(diag(real.(σ))))
    else
        σ = zeros(n_states, n_states)
        σ[state, state] = 1
    end
    return SurfaceHoppingVariables(ComponentVector(v=v, r=r, σreal=σ, σimag=zero(σ)), state)
end

function DynamicsVariables(sim::RingPolymerSimulation{<:SurfaceHopping}, v, r, state::Integer; type=:diabatic)
    n_states = sim.calculator.model.n_states
    if type == :diabatic
        Calculators.evaluate_centroid_potential!(sim.calculator, r)
        U = eigvecs(sim.calculator.centroid_potential)

        diabatic_density = zeros(n_states, n_states)
        diabatic_density[state, state] = 1
        σ = U' * diabatic_density * U
        state = sample(Weights(diag(real.(σ))))
    else
        σ = zeros(n_states, n_states)
        σ[state, state] = 1
    end
    return SurfaceHoppingVariables(ComponentVector(v=v, r=r, σreal=σ, σimag=zero(σ)), state)
end

get_velocities(u::SurfaceHoppingVariables) = get_velocities(u.x)
get_positions(u::SurfaceHoppingVariables) = get_positions(u.x)
get_quantum_subsystem(u::SurfaceHoppingVariables) = get_quantum_subsystem(u.x)
