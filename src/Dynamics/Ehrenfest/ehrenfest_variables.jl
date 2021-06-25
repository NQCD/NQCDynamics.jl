export EhrenfestVariables
export get_density_matrix

mutable struct EhrenfestVariables{T,D}  <: DynamicalVariables{T}
    x::ArrayPartition{Complex{T}, Tuple{D,D,Matrix{Complex{T}}}}
end

function EhrenfestVariables(x::ArrayPartition{T}) where {T<:AbstractFloat}
    EhrenfestVariables(ArrayPartition(x.x[1], x.x[2], Complex.(x.x[3])))
end

function EhrenfestVariables(v::AbstractArray, r::AbstractArray, n_states::Integer, state::Integer)
    σ = zeros(Complex{eltype(r)}, n_states, n_states)
    σ[state, state] = 1
    EhrenfestVariables(ArrayPartition(v, r, σ))
end

function EhrenfestVariables(sim::Simulation{<:AbstractEhrenfest}, v, r, state::Integer)
    n_states = sim.calculator.model.n_states
    Calculators.evaluate_potential!(sim.calculator, r)
    Calculators.eigen!(sim.calculator)
    U = sim.calculator.eigenvectors

    diabatic_density = zeros(Complex{eltype(r)}, n_states, n_states)
    diabatic_density[state, state] = 1
    σ = U' * diabatic_density * U

    EhrenfestVariables(ArrayPartition(v, r, σ))
end

function EhrenfestVariables(sim::RingPolymerSimulation{<:AbstractEhrenfest}, v, r, state::Integer)
    n_states = sim.calculator.model.n_states
    potential!(sim.calculator.model, sim.calculator.potential[1], dropdims(mean(r; dims=3), dims=3))
    vals, U = eigen!(sim.calculator.potential[1])

    diabatic_density = zeros(Complex{eltype(r)}, n_states, n_states)
    diabatic_density[state, state] = 1
    σ = U' * diabatic_density * U

    EhrenfestVariables(ArrayPartition(v, r, σ))
end

get_density_matrix(u::EhrenfestVariables) = u.x.x[3]
