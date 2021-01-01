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

using ....Dynamics: Method
using ....NonadiabaticMolecularDynamics

export run_monte_carlo_sampling
export MonteCarlo

"""Parameters for MonteCarlo simulation
"""
mutable struct MonteCarlo{T} <: Method
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

"""Container for storing simulation quantities
"""
mutable struct MonteCarloOutput{T<:AbstractFloat,S}
    R::Vector{S}
    acceptance::T
    energy::Vector{T}
    function MonteCarloOutput{T,S}(R0::S, length::Integer) where {T,S}
        output = [similar(R0) for _=1:length]
        energy = zeros(length)
        new(output, 0.0, energy)
    end
end
MonteCarloOutput{T}(R0::S, length::Integer) where {T,S} = MonteCarloOutput{T,S}(R0, length)

# """Constructor for Monte-Carlo system
# """
# function System{MonteCarlo}(
#     atomic_parameters::AtomicParameters{T},
#     model::Models.Model,
#     temperature::Unitful.Temperature{<:Real},
#     Δ::Dict{Symbol, T},
#     n_DoF::Integer=3;
#     passes::Real=10,
#     fix::Vector{<:Integer}=Int[]) where {T}

#     monte_carlo = MonteCarlo{T}(Δ, atomic_parameters.n_atoms, passes, fix)

#     System{MonteCarlo, T}(n_DoF, austrip(temperature),
#         atomic_parameters, model, monte_carlo)
# end

# """Constructor for the MonteCarlo RingPolymerSystem
# """
# function RingPolymerSystem{MonteCarlo}(
#     atomic_parameters::AtomicParameters{T},
#     model::Models.Model,
#     temperature::Unitful.Temperature{<:Real},
#     n_beads::Integer,
#     Δ::Dict{Symbol, T},
#     n_DoF::Integer=3;
#     passes::Real=10,
#     fix::Vector{<:Integer}=Int[],
#     quantum_nuclei::Vector{Symbol}=Symbol[]) where {T<:AbstractFloat}

#     monte_carlo = MonteCarlo{T}(Δ, atomic_parameters.n_atoms, passes, fix)

#     RingPolymerSystem{MonteCarlo, T}(n_DoF, austrip(temperature),
#         atomic_parameters, model, n_beads, quantum_nuclei, monte_carlo)
# end

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
function run_monte_carlo_sampling(sim::AbstractSimulation{<:MonteCarlo}, R0::Configuration{T}) where {T}

    Rᵢ = copy(R0) # Current positions
    Rₚ = zero(Rᵢ) # Proposed positions
    sim.method.Eᵢ = evaluate_configurational_energy(sim, Rᵢ)
    output = MonteCarloOutput{T}(Rᵢ, sim.method.steps)

    run_main_loop!(sim, Rᵢ, Rₚ, output)

    output.acceptance /= sim.method.steps
    output
end

"""
    run_main_loop!(system::System{MonteCarlo}, Rᵢ::Matrix, Rₚ::Matrix, output::MonteCarloOutput)
    
Main loop for classical systems.
"""
function run_main_loop!(sim::Simulation{<:MonteCarlo}, Rᵢ::Matrix, Rₚ::Matrix, output::MonteCarloOutput)
    @showprogress 0.1 "Sampling... " for i=1:convert(Int, sim.method.steps)
        propose_move!(sim, Rᵢ, Rₚ)
        apply_cell_boundaries!(sim, Rₚ)
        assess_proposal!(sim, Rᵢ, Rₚ, output, i)
    end
end

"""Main loop for ring polymer systems.
"""
function run_main_loop!(sim::RingPolymerSimulation{<:MonteCarlo}, Rᵢ::Array{T,3}, Rₚ::Array{T,3},
    output::MonteCarloOutput) where {T}
    @showprogress 0.1 "Sampling... " for i=1:length(sim.beads):convert(Int, sim.method.steps)
        if !all(sim.method.moveable_atoms .== 0.0)
            propose_centroid_move!(sim, Rᵢ, Rₚ)
            apply_cell_boundaries!(sim, Rₚ)
        end
        assess_proposal!(sim, Rᵢ, Rₚ, output, i)
        for j=1:length(sim.beads)-1
            propose_normal_mode_move!(sim, Rᵢ, Rₚ)
            assess_proposal!(sim, Rᵢ, Rₚ, output, i+j)
        end
    end
end

"""
    propose_move!(system::System{MonteCarlo}, Rᵢ::Matrix{T}, Rₚ::Matrix{T}) where {T<:AbstractFloat}
    
Propose simple cartesian move for a single atom.
"""
function propose_move!(sim::Simulation{MonteCarlo{T}}, Rᵢ::Matrix{T}, Rₚ::Matrix{T}) where {T<:AbstractFloat}
    Rₚ .= Rᵢ
    atom = sample(range(sim.atoms), sim.method.moveable_atoms)
    apply_random_perturbation!(sim, Rₚ, atom)
end

"""
    propose_centroid_move!(system::RingPolymerSystem{MonteCarlo}, Rᵢ::Array{T, 3}, Rₚ::Array{T, 3}) where {T<:AbstractFloat}
    
Propose a move for the ring polymer centroid for one atom.
"""
function propose_centroid_move!(sim::RingPolymerSimulation{<:MonteCarlo}, Rᵢ::Array{T, 3}, Rₚ::Array{T, 3}) where {T<:AbstractFloat}
    Rₚ .= Rᵢ
    transform_to_normal_modes!(sim.beads, Rₚ)
    atom = sample(range(sim.atoms), sim.method.moveable_atoms)
    @views apply_random_perturbation!(sim, Rₚ[:,:,1], atom) # Perturb the centroid mode
    @views if atom ∉ sim.beads.quantum_atoms # Ensure replicas are indentical if not quantum
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
function apply_random_perturbation!(sim::AbstractSimulation{<:MonteCarlo}, R::AbstractMatrix, atom::Integer)
    R[:, atom] .+= (rand(sim.DoFs) .- 0.5) .* sim.method.Δ[sim.atoms.types[atom]] / sqrt(sim.DoFs)
end

"""
    propose_normal_mode_move!(system::RingPolymerSystem{MonteCarlo}, Rᵢ::Array{T, 3}, Rₚ::Array{T, 3}) where {T}
    
Propose a move for a single normal mode for a single atom.
"""
function propose_normal_mode_move!(sim::RingPolymerSimulation{<:MonteCarlo}, Rᵢ::Array{T, 3}, Rₚ::Array{T, 3}) where {T}
    Rₚ .= Rᵢ
    transform_to_normal_modes!(sim.beads, Rₚ)
    atom = sample(sim.beads.quantum_atoms)
    mode = sample(2:length(sim.beads))
    @views sample_mode!(sim.beads, sim.atoms.masses[atom], mode, Rₚ[:, atom, mode])
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
function assess_proposal!(sim::AbstractSimulation{<:MonteCarlo}, Rᵢ, Rₚ, output, i)
    sim.method.Eₚ = evaluate_configurational_energy(sim, Rₚ)
    if acceptance_probability(sim) > rand()
        sim.method.Eᵢ = sim.method.Eₚ
        Rᵢ .= Rₚ
        output.acceptance += 1
    end
    write_output!(output, Rᵢ, sim.method.Eᵢ, i)
end

"""
    apply_cell_boundaries!(system::System{MonteCarlo}, R::Matrix)
    
Ensure the atom remains inside the simulation cell.
"""
function apply_cell_boundaries!(sim::Simulation{<:MonteCarlo}, R::Matrix)
    apply_cell_boundaries!(sim.cell, R, length(sim.atoms))
end

"""
    apply_cell_boundaries!(system::RingPolymerSystem{MonteCarlo}, R::Array{T, 3}) where {T}
    
Apply cell boundaries to the ring polymer system.

This converts to normal mode coordinates and restricts only the centroid.
This means that replicas can exit the cell but the centroid cannot.
The reasoning for this is to avoid complications in the computation of the spring potential.
Proper treatment would require accounting for the periodic cell in this function which I have
not yet done.

The modification by a factor of ``\\sqrt{N}`` is to convert to real space centroid. 
"""
function apply_cell_boundaries!(sim::RingPolymerSimulation{<:MonteCarlo}, R::Array{T, 3}) where {T}
    transform_to_normal_modes!(sim.beads, R, sim.DoFs)
    R[:,sim.beads.quantum_atoms,1] ./= sqrt(length(sim.beads))
    @views apply_cell_boundaries!(sim.cell, R[:,:,1], length(sim.atoms))
    R[:,sim.beads.quantum_atoms,1] .*= sqrt(length(sim.beads))
    transform_from_normal_modes!(sim.beads, R, sim.DoFs)
end

"""
    apply_cell_boundaries!(cell::PeriodicCell, R::Matrix, n_atoms::Integer, n_DoF::Integer)
    
Apply simple periodic boundaries 
"""
function apply_cell_boundaries!(cell::PeriodicCell, R::AbstractMatrix, n_atoms::Integer)
    @views for i=1:n_atoms
        R[:,i] .= cell.vectors \ R[:,i]
        R[:,i] .= mod1.(R[:,i], 1)
        R[:,i] .= cell.vectors * R[:,i]
    end
end
apply_cell_boundaries!(::InfiniteCell, ::AbstractMatrix, ::Integer) = nothing

"""
    acceptance_probability(system::AbstractSystem{MonteCarlo})

Return the Metropolis-Hastings acceptance probability.
"""
function acceptance_probability(sim::Simulation{<:MonteCarlo})
    acceptance_probability(sim.method.Eₚ, sim.method.Eᵢ, sim.temperature)
end
function acceptance_probability(sim::RingPolymerSimulation{<:MonteCarlo})
    acceptance_probability(sim.method.Eₚ, sim.method.Eᵢ, sim.beads.ω_n)
end
acceptance_probability(Eₚ::T, Eᵢ::T, kT::T) where {T<:AbstractFloat} = min(1, exp(-(Eₚ - Eᵢ)/kT))

"""
    write_output!(output::MonteCarloOutput, Rᵢ::Configuration, energy::AbstractFloat, i::Integer)
    
Store the current configuration and associated energy.
"""
function write_output!(output::MonteCarloOutput, Rᵢ::Configuration, energy::AbstractFloat, i::Integer)
    output.R[i] .= Rᵢ
    output.energy[i] = energy
end

end # module