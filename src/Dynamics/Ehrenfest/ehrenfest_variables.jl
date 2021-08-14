# export EhrenfestVariables
export get_quantum_subsystem
using ComponentArrays
using StructArrays

function DynamicsVariables(sim::Simulation{<:AbstractEhrenfest}, v, r, state::Integer; type=:diabatic)
    n_states = sim.calculator.model.n_states
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

function DynamicsVariables(sim::RingPolymerSimulation{<:AbstractEhrenfest}, v, r, state::Integer; type=:diabatic)
    n_states = sim.calculator.model.n_states
    if type == :diabatic
        Calculators.evaluate_centroid_potential!(sim.calculator, r)
        U = eigvecs(sim.calculator.centroid_potential)

        diabatic_density = zeros(n_states, n_states)
        diabatic_density[state, state] = 1
        σ = U' * diabatic_density * U
    else
        σ = zeros(n_states, n_states)
        σ[state, state] = 1
    end
    return ComponentVector(v=v, r=r, σreal=σ, σimag=zero(σ))
end

get_quantum_subsystem(u::ComponentVector{T}) where {T} =
    StructArray{Complex{T}}((u.σreal, u.σimag))
