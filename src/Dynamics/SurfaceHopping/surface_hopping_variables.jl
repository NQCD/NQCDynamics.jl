export SurfaceHoppingVariables
export get_density_matrix
using StatsBase: sample, Weights

mutable struct SurfaceHoppingVariables{T,D,S}  <: DynamicalVariables{T}
    x::ArrayPartition{Complex{T}, Tuple{D,D,Matrix{Complex{T}}}}
    state::S
end

function SurfaceHoppingVariables(x::ArrayPartition{T}, state) where {T<:AbstractFloat}
    SurfaceHoppingVariables(ArrayPartition(x.x[1], x.x[2], Complex.(x.x[3])), state)
end

function SurfaceHoppingVariables(v::AbstractArray, r::AbstractArray, n_states::Integer, state::Integer)
    σ = zeros(Complex{eltype(r)}, n_states, n_states)
    σ[state, state] = 1
    SurfaceHoppingVariables(ArrayPartition(v, r, σ), state)
end

function SurfaceHoppingVariables(sim::Simulation{<:SurfaceHopping}, v, r, state::Integer)
    n_states = sim.calculator.model.n_states
    Calculators.evaluate_potential!(sim.calculator, r)
    Calculators.eigen!(sim.calculator)
    U = sim.calculator.eigenvectors

    diabatic_density = zeros(Complex{eltype(r)}, n_states, n_states)
    diabatic_density[state, state] = 1
    σ = U' * diabatic_density * U
    adiabatic_state = sample(Weights(diag(real.(σ))))

    SurfaceHoppingVariables(ArrayPartition(v, r, σ), adiabatic_state)
end

get_density_matrix(u::SurfaceHoppingVariables) = u.x.x[3]
