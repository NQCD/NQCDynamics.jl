"""
    MetropolisHastings
    
Sampling of the initial conditions using the Metropolis-Hastings Markov chain Monte Carlo method.
"""
module MetropolisHastings

using StatsBase
using ProgressMeter
using Unitful
using UnitfulAtomic
using Distributions
using ....Systems
using ....Atoms
using ....Models
using ....Electronics

export run_monte_carlo_sampling

mutable struct MonteCarloOutput{T<:AbstractFloat}
    R::Vector{Matrix{T}}
    acceptance::T
    energy::Vector{T}
    function MonteCarloOutput{T}(shape, length) where {T}
        output = [Matrix{T}(undef, shape) for _=1:length]
        energy = zeros(length)
        new(output, 0.0, energy)
    end
end

function run_monte_carlo_sampling(system::System{Classical, T}, R0::Matrix{T};
    passes::Real=10, Δ::Real=1.0, fix::Vector{<:Integer}=Int[]) where {T}
    
    moveables = ones(n_atoms(system))
    moveables[fix] .= 0
    frozen = AnalyticWeights(moveables)

    Rᵢ = copy(R0)
    Rₚ = zero(Rᵢ) # Proposed positions
    energy = zeros(2)
    energy[1] = Electronics.evaluate_potential(system.model, Rᵢ)
    output = MonteCarloOutput{eltype(R0)}(size(Rᵢ), passes*n_atoms(system))

    @showprogress 0.1 "Sampling... " for i=1:convert(Int, passes*n_atoms(system))
        perform_monte_carlo_step!(system, Rᵢ, Rₚ, energy, Δ, output, i, frozen)
    end
    output.acceptance /= passes*n_atoms(system)
    output
end

function perform_monte_carlo_step!(system::System{Classical, T}, Rᵢ::Matrix{T}, Rₚ::Matrix{T},
    energy::Vector{T}, Δ::T, output::MonteCarloOutput{T}, i::Integer, frozen::AnalyticWeights) where {T<:AbstractFloat}
    propose_move!(Rᵢ, Rₚ, Δ, n_DoF(system), n_atoms(system), frozen)
    apply_cell_boundaries!(system.atomic_parameters.cell, Rₚ, n_atoms(system), n_DoF(system))
    energy[2] = Electronics.evaluate_potential(system.model, Rₚ)
    assess_proposal!(energy, Rₚ, Rᵢ, system.temperature, output)
    write_output!(output, Rᵢ, energy[1], i)
end

function propose_move!(Rᵢ::Matrix{T}, Rₚ::Matrix{T}, Δ::T, n_DoF::Integer, n_atoms::Integer,
    frozen::AnalyticWeights) where {T<:AbstractFloat}
    Rₚ .= Rᵢ
    atom = sample(1:n_atoms, frozen)
    Rₚ[:,atom] .+= (rand(n_DoF) .- 0.5) .* Δ / sqrt(n_DoF)
end

"""
    apply_cell_boundaries!(cell::PeriodicCell, R::Matrix{<:AbstractFloat}, n_atoms::Integer, n_DoF::Integer)
    
Apply simple periodic boundaries 

This works for cubic cells only.
"""
function apply_cell_boundaries!(cell::PeriodicCell, R::Matrix{<:AbstractFloat}, n_atoms::Integer, n_DoF::Integer)
    cell_vectors = austrip.(cell.vectors)
    for i=1:n_atoms
        for j=1:n_DoF # Loop ordering is very inefficient, don't need to check every time
            if cell.periodicity[j] # If periodic in this dimension
                if R[j,i] > cell_vectors[j,j]
                    R[j,i] -= cell_vectors[j,j]
                elseif R[j,i] < 0
                    R[j,i] += cell_vectors[j,j]
                end
            end
        end
    end
end
apply_cell_boundaries!(cell::InfiniteCell, R, n_atoms, n_DoF) = nothing

function assess_proposal!(energy::Vector{T}, Rₚ::Matrix{T}, Rᵢ::Matrix{T}, temperature::T, output::MonteCarloOutput{T}) where {T<:AbstractFloat}
    if acceptance_probability(energy[2], energy[1], temperature) > rand()
        energy[1] = energy[2]
        Rᵢ .= Rₚ
        output.acceptance += 1
    end
end

function acceptance_probability(energyₚ::T, energyᵢ::T, temperature::T) where {T<:AbstractFloat}
    min(1, exp(-(energyₚ - energyᵢ)/temperature))
end

function write_output!(output::MonteCarloOutput, Rᵢ::Matrix{T}, energy::T, i::Integer) where {T<:AbstractFloat}
    output.R[i] .= Rᵢ
    output.energy[i] = energy
end

end # module