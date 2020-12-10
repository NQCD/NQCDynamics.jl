"""
    MetropolisHastings
    
Sampling of the initial conditions using the Metropolis-Hastings Markov chain Monte Carlo method.
"""
module MetropolisHastings

using ProgressMeter
using Unitful
using UnitfulAtomic
using Distributions
using ....Systems
using ....Atoms
using ....Models

export run_monte_carlo_sampling

mutable struct MonteCarloOutput{T<:AbstractFloat}
    R::Vector{Matrix{T}}
    acceptance::T
    function MonteCarloOutput{T}(shape, length) where {T}
        output = [Matrix{T}(undef, shape) for _=1:length]
        new(output, 0.0)
    end
end

function run_monte_carlo_sampling(system::System{Classical, T}, R0::Matrix{T}; passes::Real=10, Δ::Real=1.0) where {T}

    Rᵢ = copy(R0)
    Rₚ = zero(Rᵢ) # Proposed positions
    energy = zeros(2)
    energy[1] = system.model.get_V0(Rᵢ)
    output = MonteCarloOutput{eltype(R0)}(size(Rᵢ), passes*n_atoms(system))

    @showprogress 0.1 "Sampling... " for i=1:convert(Int, passes*n_atoms(system))
        perform_monte_carlo_step!(system, Rᵢ, Rₚ, energy, Δ, output, i)
    end
    output.acceptance /= passes
    output
end

function perform_monte_carlo_step!(system::System{Classical, T}, Rᵢ::Matrix{T}, Rₚ::Matrix{T},
    energy::Vector{T}, Δ::T, output::MonteCarloOutput{T}, i::Integer) where {T<:AbstractFloat}
    propose_move!(Rᵢ, Rₚ, Δ, n_DoF(system), n_atoms(system))
    energy[2] = system.model.get_V0(Rₚ)
    assess_proposal!(energy, Rₚ, Rᵢ, system.temperature, output)
    write_output!(output.R[i], Rᵢ)
end

function propose_move!(Rᵢ::Matrix{T}, Rₚ::Matrix{T}, Δ::T, n_DoF::Integer, n_atoms::Integer) where {T<:AbstractFloat}
    Rₚ .= Rᵢ
    atom = sample(1:n_atoms)
    Rₚ[:,atom] .+= (rand(n_DoF) .- 0.5) .* Δ / sqrt(n_DoF)
end

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

write_output!(output, Rₚ) = output .= Rₚ

end # module