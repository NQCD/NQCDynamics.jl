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
using DocStringExtensions

using ....NonadiabaticMolecularDynamics

export run_monte_carlo_sampling
export MonteCarlo
export PathIntegralMonteCarlo

abstract type MonteCarloParameters end

"""
$(TYPEDEF)

Parameters for Monte carlo simulations.
"""
mutable struct MonteCarlo{T} <: MonteCarloParameters
    Eᵢ::T
    Eₚ::T
    Δ::Dict{Symbol, T}
    steps::UInt
    moveable_atoms::AnalyticWeights
    extra_function::Function
    function MonteCarlo{T}(Δ::Dict{Symbol, T},
        n_atoms::Real, passes::Real, fix::Vector{<:Integer}, extra::Function) where {T<:AbstractFloat}

        steps = n_atoms * passes

        moveables = ones(n_atoms)
        moveables[fix] .= 0
        moveable_atoms = AnalyticWeights(moveables)

        new(0.0, 0.0, Δ, steps, moveable_atoms, extra)
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
    acceptance::Dict{Symbol, T}
    total_moves::Dict{Symbol, UInt}
    energy::Vector{T}
    function MonteCarloOutput(R0::S, atoms::Atoms{N,T}) where {N,T,S}
        output = typeof(R0)[]
        acceptance = Dict{Symbol, T}()
        total_moves = Dict{Symbol, UInt}()
        for element in atoms.types
            acceptance[element] = 0.0
            total_moves[element] = 0
        end
        energy = T[]
        new{T,S}(output, acceptance, total_moves, energy)
    end
end

"""
    run_monte_carlo_sampling(sim::AbstractSimulation, monte::MonteCarloParameters, R0)
    
Perform Monte Carlo sampling for the system defined by the `sim` and `monte` parameters.

From the initial positions specified `R0` the system will be explored using the
Metropolis-Hastings algorithm.
"""
function run_monte_carlo_sampling(sim::Simulation, R0::Matrix{T},
    Δ::Dict{Symbol,T}, passes::Real; fix::Vector{<:Integer}=Int[], extra_function::Function=x->true) where {T}

    Rᵢ = copy(R0) # Current positions
    Rₚ = zero(Rᵢ) # Proposed positions
    monte = MonteCarlo{T}(Δ, length(sim.atoms), passes, fix, extra_function)
    monte.Eᵢ = evaluate_potential_energy(sim, Rᵢ)
    output = MonteCarloOutput(Rᵢ, sim.atoms)

    run_main_loop!(sim, monte, Rᵢ, Rₚ, output)

    for key in keys(output.acceptance)
        output.acceptance[key] /= output.total_moves[key]
    end
    output
end

function run_monte_carlo_sampling(sim::RingPolymerSimulation, R0::Array{T,3},
    Δ::Dict{Symbol,T}, passes::Real; fix::Vector{<:Integer}=Int[], segment::Real=1) where {T}

    Rᵢ = copy(R0) # Current positions
    Rₚ = zero(Rᵢ) # Proposed positions
    monte = PathIntegralMonteCarlo{T}(Δ, length(sim.atoms), passes, fix, segment, length(sim.beads))
    monte.Eᵢ = evaluate_potential_energy(sim, Rᵢ)
    output = MonteCarloOutput(Rᵢ, sim.atoms)

    run_main_loop!(sim, monte, Rᵢ, Rₚ, output)

    for key in keys(output.acceptance)
        output.acceptance[key] /= output.total_moves[key]
    end
    output
end

"""
    run_main_loop!(system::System{MonteCarlo}, Rᵢ::Matrix, Rₚ::Matrix, output::MonteCarloOutput)
    
Main loop for classical systems.
"""
function run_main_loop!(sim::Simulation, monte::MonteCarlo, Rᵢ::Matrix, Rₚ::Matrix, output::MonteCarloOutput)
    @showprogress 0.1 "Sampling... " for i=1:convert(Int, monte.steps)
        atom = sample(range(sim.atoms), monte.moveable_atoms)
        propose_move!(sim, monte, Rᵢ, Rₚ, atom)
        @views for i in axes(Rₚ, 2) # atoms
            if monte.moveable_atoms[i] > 0
                apply_cell_boundaries!(sim.cell, Rₚ[:,i])
            end
        end
        assess_proposal!(sim, monte, Rᵢ, Rₚ, output, atom)
    end
end

"""Main loop for ring polymer systems.
"""
function run_main_loop!(sim::RingPolymerSimulation, monte::PathIntegralMonteCarlo, Rᵢ::Array{T,3}, Rₚ::Array{T,3},
    output::MonteCarloOutput) where {T}
    @showprogress 0.1 "Sampling... " for i=1:convert(Int, monte.steps)
        if !all(monte.moveable_atoms .== 0.0)
            atom = sample(range(sim.atoms), monte.moveable_atoms)
            propose_centroid_move!(sim, monte, Rᵢ, Rₚ, atom)
            apply_cell_boundaries!(sim.cell, Rₚ, sim.beads)
            assess_proposal!(sim, monte, Rᵢ, Rₚ, output, atom)
            if atom ∈ sim.beads.quantum_atoms
                for j=1:monte.pass_length
                    propose_normal_mode_move!(sim, monte, Rᵢ, Rₚ, atom)
                    assess_proposal!(sim, monte, Rᵢ, Rₚ, output, atom)
                end
            end
        else
            atom = sample(sim.beads.quantum_atoms)
            for j=1:monte.pass_length
                propose_normal_mode_move!(sim, monte, Rᵢ, Rₚ, atom)
                assess_proposal!(sim, monte, Rᵢ, Rₚ, output, atom)
            end
        end
        assess_proposal!(sim, monte, Rᵢ, Rₚ, output, atom)
    end
end

"""
    propose_move!(system::System{MonteCarlo}, Rᵢ::Matrix{T}, Rₚ::Matrix{T}) where {T<:AbstractFloat}
    
Propose simple cartesian move for a single atom.
"""
function propose_move!(sim::Simulation, monte::MonteCarlo, Rᵢ::Matrix, Rₚ::Matrix, atom::Integer)
    Rₚ .= Rᵢ
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
function assess_proposal!(sim::AbstractSimulation, monte::MonteCarloParameters, Rᵢ, Rₚ, output, atom::Integer)
    monte.Eₚ = evaluate_potential_energy(sim, Rₚ)
    output.total_moves[sim.atoms.types[atom]] += 1
    if monte.extra_function(Rₚ)
        if acceptance_probability(sim, monte) > rand()
            monte.Eᵢ = monte.Eₚ
            Rᵢ .= Rₚ
            output.acceptance[sim.atoms.types[atom]] += 1
        end
    end
    write_output!(output, Rᵢ, monte.Eᵢ)
end

"""
    acceptance_probability(system::AbstractSystem{MonteCarlo})

Return the Metropolis-Hastings acceptance probability.
"""
function acceptance_probability(sim::Simulation, monte::MonteCarlo)
    acceptance_probability(monte.Eₚ, monte.Eᵢ, get_temperature(sim))
end
function acceptance_probability(sim::RingPolymerSimulation, monte::PathIntegralMonteCarlo)
    acceptance_probability(monte.Eₚ, monte.Eᵢ, sim.beads.ω_n)
end
acceptance_probability(Eₚ::T, Eᵢ::T, kT::T) where {T<:AbstractFloat} = min(1, exp(-(Eₚ - Eᵢ)/kT))

"""
    write_output!(output::MonteCarloOutput, Rᵢ::AbstractArray, energy::AbstractFloat, i::Integer)
    
Store the current configuration and associated energy.
"""
function write_output!(output::MonteCarloOutput, Rᵢ::AbstractArray, energy::AbstractFloat)
    push!(output.R, copy(Rᵢ))
    push!(output.energy, copy(energy))
end

end # module