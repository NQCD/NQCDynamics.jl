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
export MonteCarlo

"""Parameters for MonteCarlo simulation"""
mutable struct MonteCarlo{T<:AbstractFloat} <: Systems.DynamicsParameters
    Eᵢ::T
    Eₚ::T
    Δ::Dict{Symbol, T}
    steps::UInt
    moveable_atoms::AnalyticWeights
end

"""Container for storing simulation quantities"""
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

"""Constructor for Monte-Carlo system"""
function System{MonteCarlo}(
    atomic_parameters::AtomicParameters{T},
    model::Models.Model,
    temperature::Unitful.Temperature{<:Real},
    Δ::Dict{Symbol, T},
    n_DoF::Integer=3;
    passes::Real=10,
    fix::Vector{<:Integer}=Int[]) where {T}

    moveables = ones(atomic_parameters.n_atoms)
    moveables[fix] .= 0
    moveable_atoms = AnalyticWeights(moveables)

    steps = passes * atomic_parameters.n_atoms
    monte_carlo = MonteCarlo{T}(0.0, 0.0, Δ, steps, moveable_atoms)

    System{MonteCarlo, T}(n_DoF, austrip(temperature),
        atomic_parameters, model, monte_carlo)
end

"""
    run_monte_carlo_sampling(system::System{MonteCarlo}, R0)
    
Run the simulation defined by the MonteCarlo system.
"""
function run_monte_carlo_sampling(system::System{MonteCarlo}, R0)

    Rᵢ = copy(R0) # Current positions
    Rₚ = zero(Rᵢ) # Proposed positions
    system.dynamics.Eᵢ = Electronics.evaluate_potential(system.model, Rᵢ)
    output = MonteCarloOutput{eltype(R0)}(size(Rᵢ), system.dynamics.steps)

    @showprogress 0.1 "Sampling... " for i=1:convert(Int, system.dynamics.steps)
        perform_monte_carlo_step!(system, Rᵢ, Rₚ, output, i)
    end

    output.acceptance /= system.dynamics.steps
    output
end

"""
    perform_monte_carlo_step!(system::System{MonteCarlo}, Rᵢ::Matrix{T}, Rₚ::Matrix{T},
    
Perform a single step of the MetropolisHastings algorithm.
"""
function perform_monte_carlo_step!(system::System{MonteCarlo}, Rᵢ::Matrix{T}, Rₚ::Matrix{T},
    output::MonteCarloOutput{T}, i::Integer) where {T<:AbstractFloat}
    propose_move!(system, Rᵢ, Rₚ)
    apply_cell_boundaries!(system, Rₚ)
    system.dynamics.Eₚ = Electronics.evaluate_potential(system.model, Rₚ)
    assess_proposal!(system, Rₚ, Rᵢ, output)
    write_output!(output, Rᵢ, system.dynamics.Eᵢ, i)
end

function propose_move!(system::System{MonteCarlo}, Rᵢ::Matrix{T}, Rₚ::Matrix{T}) where {T<:AbstractFloat}
    Rₚ .= Rᵢ
    atom = sample(1:n_atoms(system), system.dynamics.moveable_atoms)
    Rₚ[:,atom] .+= (rand(n_DoF(system)) .- 0.5) .* system.dynamics.Δ[system.atomic_parameters.atom_types[atom]] / sqrt(n_DoF(system))
end

function assess_proposal!(system::System{MonteCarlo}, Rₚ::Matrix{T}, Rᵢ::Matrix{T}, output::MonteCarloOutput{T}) where {T<:AbstractFloat}
    if acceptance_probability(system) > rand()
        system.dynamics.Eᵢ = system.dynamics.Eₚ
        Rᵢ .= Rₚ
        output.acceptance += 1
    end
end

function apply_cell_boundaries!(system::System{MonteCarlo}, R::Matrix{<:AbstractFloat})
    apply_cell_boundaries!(system.atomic_parameters.cell, R, n_atoms(system), n_DoF(system))
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

function acceptance_probability(system::System{MonteCarlo})
    acceptance(system.dynamics.Eₚ, system.dynamics.Eᵢ, system.temperature)
end
acceptance(Eₚ::T, Eᵢ::T, kT::T) where {T<:AbstractFloat} = min(1, exp(-(Eₚ - Eᵢ)/kT))

function write_output!(output::MonteCarloOutput, Rᵢ::Matrix{T}, energy::T, i::Integer) where {T<:AbstractFloat}
    output.R[i] .= Rᵢ
    output.energy[i] = energy
end

end # module