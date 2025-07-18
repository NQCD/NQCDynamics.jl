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
using StatsBase: AnalyticWeights, sample
using ProgressMeter: @showprogress
using Distributions: Normal
using RingPolymerArrays: RingPolymerArrays

using NQCBase: NQCBase, Atoms
using NQCCalculators
using NQCDynamics:
    NQCDynamics,
    AbstractSimulation,
    Simulation,
    RingPolymerSimulation,
    RingPolymers,
    Estimators,
    DynamicsUtils,
    ndofs
export run_monte_carlo_sampling

abstract type MonteCarloParameters end

"""
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
    extra_function::Function
    function PathIntegralMonteCarlo{T}(Δ::Dict{Symbol, T},
        n_atoms::Real, passes::Real, fix::Vector{<:Integer},
        segment::Real, n_beads::Integer, extra::Function) where {T<:AbstractFloat}

        steps = n_atoms * passes

        moveables = ones(n_atoms)
        moveables[fix] .= 0
        moveable_atoms = AnalyticWeights(moveables)

        segment_length = max(1, round(Int, segment * n_beads))
        pass_length = n_beads ÷ segment_length

        new(0.0, 0.0, Δ, steps, moveable_atoms, segment_length, pass_length, extra)
    end
end

"""Container for storing simulation quantities
"""
mutable struct MonteCarloOutput{T<:AbstractFloat,S}
    R::Vector{S}
    acceptance::Dict{Symbol, T}
    total_moves::Dict{Symbol, UInt}
    energy::Vector{T}
    function MonteCarloOutput(R0::S, atoms::Atoms{T}) where {T,S}
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
    NQCCalculators.update_cache!(sim.cache, Rᵢ)
    monte.Eᵢ = DynamicsUtils.classical_potential_energy(sim, Rᵢ)
    output = MonteCarloOutput(Rᵢ, sim.atoms)

    run_main_loop!(sim, monte, Rᵢ, Rₚ, output)

    for key in keys(output.acceptance)
        output.acceptance[key] /= output.total_moves[key]
    end
    NQCCalculators.update_cache!(sim.cache, R0)
    return output
end

function run_monte_carlo_sampling(sim::RingPolymerSimulation, R0::AbstractArray{T,3},
    Δ::Dict{Symbol,T}, passes::Real;
    fix::Vector{<:Integer}=Int[], segment::Real=1, extra_function::Function=x->true) where {T}

    Rᵢ = copy(R0) # Current positions
    Rₚ = zero(Rᵢ) # Proposed positions
    monte = PathIntegralMonteCarlo{T}(Δ, length(sim.atoms), passes, fix, segment, length(sim.beads), extra_function)
    NQCCalculators.update_cache!(sim.cache, Rᵢ)
    monte.Eᵢ = DynamicsUtils.classical_potential_energy(sim, Rᵢ)
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
                NQCBase.apply_cell_boundaries!(sim.cell, Rₚ[:,i])
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
            NQCBase.apply_cell_boundaries!(sim.cell, Rₚ, sim.beads)
            assess_proposal!(sim, monte, Rᵢ, Rₚ, output, atom)
            if length(sim.beads) > 1
                if atom ∈ sim.beads.quantum_atoms
                    for j=1:monte.pass_length
                        propose_normal_mode_move!(sim, monte, Rᵢ, Rₚ, atom)
                        assess_proposal!(sim, monte, Rᵢ, Rₚ, output, atom)
                    end
                end
            end
        else
            if length(sim.beads) > 1
                atom = sample(sim.beads.quantum_atoms)
                for j=1:monte.pass_length
                    propose_normal_mode_move!(sim, monte, Rᵢ, Rₚ, atom)
                    assess_proposal!(sim, monte, Rᵢ, Rₚ, output, atom)
                end
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
    apply_random_perturbation!(sim.atoms, monte, Rₚ, atom, ndofs(sim))
end

"""
    propose_centroid_move!(system::RingPolymerSystem{MonteCarlo}, Rᵢ::Array{T, 3}, Rₚ::Array{T, 3}) where {T<:AbstractFloat}
    
Propose a move for the ring polymer centroid for one atom.
"""
function propose_centroid_move!(sim::RingPolymerSimulation, monte::PathIntegralMonteCarlo, Rᵢ::Array{T, 3}, Rₚ::Array{T, 3}, atom::Integer) where {T<:AbstractFloat}
    Rₚ .= Rᵢ
    RingPolymerArrays.transform_to_normal_modes!(Rₚ, sim.beads.transformation)
    @views apply_random_perturbation!(sim.atoms, monte, Rₚ[:,:,1], atom, ndofs(sim)) # Perturb the centroid mode
    @views if atom ∉ sim.beads.quantum_atoms # Ensure replicas are identical if not quantum
        for i=2:length(sim.beads)
            Rₚ[:, atom, i] .= Rₚ[:, atom, 1]
        end
    end
    RingPolymerArrays.transform_from_normal_modes!(Rₚ, sim.beads.transformation)
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
    RingPolymerArrays.transform_to_normal_modes!(Rₚ, sim.beads.transformation)
    for _=1:monte.segment_length
        mode = sample(2:length(sim.beads))
        @views sample_mode!(sim.beads, sim.atoms.masses[i], mode, Rₚ[:, i, mode])
    end
    RingPolymerArrays.transform_from_normal_modes!(Rₚ, sim.beads.transformation)
end

function sample_mode!(beads::RingPolymers.RingPolymerParameters, mass::AbstractFloat, mode::Integer, R::AbstractArray)
    σ = sqrt(beads.ω_n / mass) / (beads.normal_mode_springs[mode] * mass)
    rand!(Normal(0, σ), R)
end

"""
    assess_proposal!(system::AbstractSystem{MonteCarlo}, Rᵢ, Rₚ, output, i)
    
Update the energy, check for acceptance, and update the output. 
"""
function assess_proposal!(sim::AbstractSimulation, monte::MonteCarloParameters, Rᵢ, Rₚ, output, atom::Integer)
    NQCCalculators.update_cache!(sim.cache, Rₚ)
    monte.Eₚ = DynamicsUtils.classical_potential_energy(sim, Rₚ)
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
    acceptance_probability(monte.Eₚ, monte.Eᵢ, NQCDynamics.get_temperature(sim))
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
