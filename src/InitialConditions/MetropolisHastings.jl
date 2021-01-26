"""
    MetropolisHastings
    
Sampling of the initial conditions using the Metropolis-Hastings Markov chain Monte Carlo method.

Included within is the ability to sample the canonical distribution for adiabatic classical and
ring polymer systems.

Usage involves creating an instance of an `AbstractSystem{MonteCarlo}` and calling
`run_monte_carlo_sampling`.
"""
module MetropolisHastings

using Random
using StatsBase
using ProgressMeter
using Unitful
using UnitfulAtomic
using Distributions

using ....NonadiabaticMolecularDynamics

export run_monte_carlo_sampling
export MonteCarlo
export PathIntegralMonteCarlo

abstract type MonteCarloParameters end

"""Parameters for MonteCarlo simulation
"""
mutable struct MonteCarlo{T} <: MonteCarloParameters
    Eᵢ::T
    Eₚ::T
    Δ::Dict{Symbol, T}
    steps::UInt
    moveable_atoms::AnalyticWeights
    function MonteCarlo{T}(Δ::Dict{Symbol, T},
        n_atoms::Real, passes::Real, fix::Vector{<:Integer}) where {T<:AbstractFloat}

        steps = n_atoms * passes

        moveables = ones(n_atoms)
        moveables[fix] .= 0
        moveable_atoms = AnalyticWeights(moveables)

        new(0.0, 0.0, Δ, steps, moveable_atoms)
    end
end

mutable struct PathIntegralMonteCarlo{T} <: MonteCarloParameters
    Eᵢ::T
    Eₚ::T
    Δ::Dict{Symbol, T}
    steps::UInt
    moveable_atoms::AnalyticWeights
    segment_length::UInt
    pass_length::UInt
    function PathIntegralMonteCarlo{T}(Δ::Dict{Symbol, T},
        n_atoms::Real, passes::Real, fix::Vector{<:Integer},
        segment::Real, n_beads::Integer) where {T<:AbstractFloat}

        steps = n_atoms * passes

        moveables = ones(n_atoms)
        moveables[fix] .= 0
        moveable_atoms = AnalyticWeights(moveables)

        segment_length = max(1, round(Int, segment * n_beads))
        pass_length = n_beads ÷ segment_length

        new(0.0, 0.0, Δ, steps, moveable_atoms, segment_length, pass_length)
    end
end

"""Container for storing simulation quantities
"""
mutable struct MonteCarloOutput{T<:AbstractFloat,S}
    R::Vector{S}
    acceptance::T
    energy::Vector{T}
    function MonteCarloOutput{T,S}(R0::S) where {T,S}
        output = typeof(R0)[]
        energy = T[]
        new(output, 0.0, energy)
    end
end
function MonteCarloOutput{T}(R0::S) where {T,S}
    MonteCarloOutput{T,S}(R0)
end

"""
    Configuration{T} = Union{Matrix{T}, Array{T, 3}}

Type union for the two allowed positions types.

A Matrix for a simple classical system and a Array{T, 3} for ring polymer systems.
"""
const Configuration{T} = Union{Matrix{T}, Array{T, 3}}

"""
    function run_monte_carlo_sampling(system::AbstractSystem{MonteCarlo}, R0::Configuration{T}) where {T}
    
Run the simulation defined by the MonteCarlo system.
"""
function run_monte_carlo_sampling(sim::AbstractSimulation, monte::MonteCarloParameters, R0::Configuration{T}) where {T}

    Rᵢ = copy(R0) # Current positions
    Rₚ = zero(Rᵢ) # Proposed positions
    monte.Eᵢ = evaluate_configurational_energy(sim, Rᵢ)
    output = MonteCarloOutput{T}(Rᵢ)

    run_main_loop!(sim, monte, Rᵢ, Rₚ, output)

    output.acceptance /= length(output.R)
    output
end

"""
    run_main_loop!(system::System{MonteCarlo}, Rᵢ::Matrix, Rₚ::Matrix, output::MonteCarloOutput)
    
Main loop for classical systems.
"""
function run_main_loop!(sim::Simulation, monte::MonteCarlo, Rᵢ::Matrix, Rₚ::Matrix, output::MonteCarloOutput)
    @showprogress 0.1 "Sampling... " for i=1:convert(Int, monte.steps)
        propose_move!(sim, monte, Rᵢ, Rₚ)
        apply_cell_boundaries!(sim.cell, Rₚ)
        assess_proposal!(sim, monte, Rᵢ, Rₚ, output)
    end
end

"""Main loop for ring polymer systems.
"""
function run_main_loop!(sim::RingPolymerSimulation, monte::PathIntegralMonteCarlo, Rᵢ::Array{T,3}, Rₚ::Array{T,3},
    output::MonteCarloOutput) where {T}
    @showprogress 0.1 "Sampling... " for i=1:monte.pass_length+1:convert(Int, monte.steps*(monte.pass_length+1))
        if !all(monte.moveable_atoms .== 0.0)
            atom = sample(range(sim.atoms), monte.moveable_atoms)
            propose_centroid_move!(sim, monte, Rᵢ, Rₚ, atom)
            apply_cell_boundaries!(sim.cell, Rₚ, sim.beads)
            assess_proposal!(sim, monte, Rᵢ, Rₚ, output)
            if atom ∈ sim.beads.quantum_atoms
                for j=1:monte.pass_length
                    propose_normal_mode_move!(sim, monte, Rᵢ, Rₚ, atom)
                    assess_proposal!(sim, monte, Rᵢ, Rₚ, output)
                end
            end
        else
            atom = sample(sim.beads.quantum_atoms)
            for j=1:monte.pass_length
                propose_normal_mode_move!(sim, monte, Rᵢ, Rₚ, atom)
                assess_proposal!(sim, monte, Rᵢ, Rₚ, output)
            end
        end
        assess_proposal!(sim, monte, Rᵢ, Rₚ, output)
    end
end

"""
    propose_move!(system::System{MonteCarlo}, Rᵢ::Matrix{T}, Rₚ::Matrix{T}) where {T<:AbstractFloat}
    
Propose simple cartesian move for a single atom.
"""
function propose_move!(sim::Simulation, monte::MonteCarlo, Rᵢ::Matrix, Rₚ::Matrix)
    Rₚ .= Rᵢ
    atom = sample(range(sim.atoms), monte.moveable_atoms)
    apply_random_perturbation!(sim.atoms, monte, Rₚ, atom, sim.DoFs)
end

"""
    propose_centroid_move!(system::RingPolymerSystem{MonteCarlo}, Rᵢ::Array{T, 3}, Rₚ::Array{T, 3}) where {T<:AbstractFloat}
    
Propose a move for the ring polymer centroid for one atom.
"""
function propose_centroid_move!(sim::RingPolymerSimulation, monte::PathIntegralMonteCarlo, Rᵢ::Array{T, 3}, Rₚ::Array{T, 3}, atom::Integer) where {T<:AbstractFloat}
    Rₚ .= Rᵢ
    transform_to_normal_modes!(sim.beads, Rₚ)
    @views apply_random_perturbation!(sim.atoms, monte, Rₚ[:,:,1], atom, sim.DoFs) # Perturb the centroid mode
    @views if atom ∉ sim.beads.quantum_atoms # Ensure replicas are identical if not quantum
        for i=2:length(sim.beads)
            Rₚ[:, atom, i] .= Rₚ[:, atom, 1]
        end
    end
    transform_from_normal_modes!(sim.beads, Rₚ)
end

"""
    apply_random_perturbation!(system::AbstractSystem{MonteCarlo}, R::AbstractMatrix, atom::Integer)
    
Randomly perturb the xyz coordinates of a single atom.
"""
function apply_random_perturbation!(atoms::Atoms, monte::MonteCarloParameters, R::AbstractMatrix, atom::Integer, DoFs::Integer)
    R[:, atom] .+= (rand(DoFs) .- 0.5) .* monte.Δ[atoms.types[atom]] / sqrt(DoFs)
end

"""
    propose_normal_mode_move!(system::RingPolymerSystem{MonteCarlo}, Rᵢ::Array{T, 3}, Rₚ::Array{T, 3}) where {T}
    
Propose a move for a single normal mode for a single atom.
"""
function propose_normal_mode_move!(sim::RingPolymerSimulation, monte::PathIntegralMonteCarlo, Rᵢ::Array{T,3},  Rₚ::Array{T, 3}, i::Integer) where {T}
    Rₚ .= Rᵢ
    transform_to_normal_modes!(sim.beads, Rₚ)
    for _=1:monte.segment_length
        mode = sample(2:length(sim.beads))
        @views sample_mode!(sim.beads, sim.atoms.masses[i], mode, Rₚ[:, i, mode])
    end
    transform_from_normal_modes!(sim.beads, Rₚ)
end

function sample_mode!(beads::RingPolymerParameters, mass::AbstractFloat, mode::Integer, R::AbstractArray)
    σ = sqrt(beads.ω_n / mass) / (beads.normal_mode_springs[mode] * mass)
    rand!(Normal(0, σ), R)
end

"""
    assess_proposal!(system::AbstractSystem{MonteCarlo}, Rᵢ, Rₚ, output, i)
    
Update the energy, check for acceptance, and update the output. 
"""
function assess_proposal!(sim::AbstractSimulation, monte::MonteCarloParameters, Rᵢ, Rₚ, output)
    monte.Eₚ = evaluate_configurational_energy(sim, Rₚ)
    if acceptance_probability(sim, monte) > rand()
        monte.Eᵢ = monte.Eₚ
        Rᵢ .= Rₚ
        output.acceptance += 1
    end
    write_output!(output, Rᵢ, monte.Eᵢ)
end

"""
    acceptance_probability(system::AbstractSystem{MonteCarlo})

Return the Metropolis-Hastings acceptance probability.
"""
function acceptance_probability(sim::Simulation, monte::MonteCarlo)
    acceptance_probability(monte.Eₚ, monte.Eᵢ, sim.temperature)
end
function acceptance_probability(sim::RingPolymerSimulation, monte::PathIntegralMonteCarlo)
    acceptance_probability(monte.Eₚ, monte.Eᵢ, sim.beads.ω_n)
end
acceptance_probability(Eₚ::T, Eᵢ::T, kT::T) where {T<:AbstractFloat} = min(1, exp(-(Eₚ - Eᵢ)/kT))

"""
    write_output!(output::MonteCarloOutput, Rᵢ::Configuration, energy::AbstractFloat, i::Integer)
    
Store the current configuration and associated energy.
"""
function write_output!(output::MonteCarloOutput, Rᵢ::Configuration, energy::AbstractFloat)
    # output.R[i] .= Rᵢ
    push!(output.R, copy(Rᵢ))
    push!(output.energy, copy(energy))
    # output.energy[i] = energy
end

end # module