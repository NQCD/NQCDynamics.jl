export EhrenfestVariables
export get_quantum_subsystem

const EhrenfestVariables{T,D} = ArrayPartition{Complex{T}, Tuple{D,D,Matrix{Complex{T}}}}

function EhrenfestVariables(v::AbstractArray, r::AbstractArray, n_states::Integer, state::Integer)
    σ = zeros(Complex{eltype(r)}, n_states, n_states)
    σ[state, state] = 1
    ArrayPartition(v, r, σ)
end

function EhrenfestVariables(sim::Simulation{<:AbstractEhrenfest}, v, r, state::Integer)
    n_states = sim.calculator.model.n_states
    Calculators.evaluate_potential!(sim.calculator, r)
    Calculators.eigen!(sim.calculator)
    U = sim.calculator.eigenvectors

    diabatic_density = zeros(Complex{eltype(r)}, n_states, n_states)
    diabatic_density[state, state] = 1
    σ = U' * diabatic_density * U

    ArrayPartition(v, r, σ)
end

function EhrenfestVariables(sim::RingPolymerSimulation{<:AbstractEhrenfest}, v, r, state::Integer)
    n_states = sim.calculator.model.n_states
    Calculators.evaluate_centroid_potential!(sim.calculator, r)
    U = eigvecs(sim.calculator.centroid_potential)

    diabatic_density = zeros(Complex{eltype(r)}, n_states, n_states)
    diabatic_density[state, state] = 1
    σ = U' * diabatic_density * U

    ArrayPartition(v, r, σ)
end

get_quantum_subsystem(u::EhrenfestVariables) = u.x[3]
