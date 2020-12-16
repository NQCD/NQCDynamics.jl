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
using ....Systems
using ....Atoms
using ....Models
using ....Electronics
using ....Dynamics

export run_monte_carlo_sampling
export MonteCarlo

"""Parameters for MonteCarlo simulation
"""
mutable struct MonteCarlo{T<:AbstractFloat} <: Systems.DynamicsParameters
    Eᵢ::T
    Eₚ::T
    Δ::Dict{Symbol, T}
    steps::UInt
    moveable_atoms::AnalyticWeights
    function MonteCarlo{T}(Δ::Dict{Symbol, T},
        n_atoms::Real, passes::Real, fix::Vector{<:Integer}) where {T}

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

"""Constructor for Monte-Carlo system
"""
function System{MonteCarlo}(
    atomic_parameters::AtomicParameters{T},
    model::Models.Model,
    temperature::Unitful.Temperature{<:Real},
    Δ::Dict{Symbol, T},
    n_DoF::Integer=3;
    passes::Real=10,
    fix::Vector{<:Integer}=Int[]) where {T}

    monte_carlo = MonteCarlo{T}(Δ, atomic_parameters.n_atoms, passes, fix)

    System{MonteCarlo, T}(n_DoF, austrip(temperature),
        atomic_parameters, model, monte_carlo)
end

"""Constructor for the MonteCarlo RingPolymerSystem
"""
function RingPolymerSystem{MonteCarlo}(
    atomic_parameters::AtomicParameters{T},
    model::Models.Model,
    temperature::Unitful.Temperature{<:Real},
    n_beads::Integer,
    Δ::Dict{Symbol, T},
    n_DoF::Integer=3;
    passes::Real=10,
    fix::Vector{<:Integer}=Int[],
    quantum_nuclei::Vector{Symbol}=Symbol[]) where {T<:AbstractFloat}

    monte_carlo = MonteCarlo{T}(Δ, atomic_parameters.n_atoms, passes, fix)

    RingPolymerSystem{MonteCarlo, T}(n_DoF, austrip(temperature),
        atomic_parameters, model, n_beads, quantum_nuclei, monte_carlo)
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
function run_monte_carlo_sampling(system::AbstractSystem{MonteCarlo}, R0::Configuration{T}) where {T}

    Rᵢ = copy(R0) # Current positions
    Rₚ = zero(Rᵢ) # Proposed positions
    system.dynamics.Eᵢ = evaluate_configurational_energy(system, Rᵢ)
    output = MonteCarloOutput{T}(Rᵢ, system.dynamics.steps)

    run_main_loop!(system, Rᵢ, Rₚ, output)

    output.acceptance /= system.dynamics.steps
    output
end

"""
    run_main_loop!(system::System{MonteCarlo}, Rᵢ::Matrix, Rₚ::Matrix, output::MonteCarloOutput)
    
Main loop for classical systems.
"""
function run_main_loop!(system::System{MonteCarlo}, Rᵢ::Matrix, Rₚ::Matrix, output::MonteCarloOutput)
    @showprogress 0.1 "Sampling... " for i=1:convert(Int, system.dynamics.steps)
        propose_move!(system, Rᵢ, Rₚ)
        apply_cell_boundaries!(system, Rₚ)
        assess_proposal!(system, Rₚ, Rᵢ, output, i)
    end
end

"""Main loop for ring polymer systems.
"""
function run_main_loop!(system::RingPolymerSystem{MonteCarlo}, Rᵢ::Array{T,3}, Rₚ::Array{T,3},
    output::MonteCarloOutput) where {T}
    @showprogress 0.1 "Sampling... " for i=1:n_beads(system):convert(Int, system.dynamics.steps)
        propose_centroid_move!(system, Rᵢ, Rₚ)
        apply_cell_boundaries!(system, Rₚ)
        assess_proposal!(system, Rᵢ, Rₚ, output, i)
        for j=1:n_beads(system)-1
            propose_normal_mode_move!(system, Rᵢ, Rₚ)
            assess_proposal!(system, Rᵢ, Rₚ, output, i+j)
        end
    end
end

"""
    propose_move!(system::System{MonteCarlo}, Rᵢ::Matrix{T}, Rₚ::Matrix{T}) where {T<:AbstractFloat}
    
Propose simple cartesian move for a single atom.
"""
function propose_move!(system::System{MonteCarlo}, Rᵢ::Matrix{T}, Rₚ::Matrix{T}) where {T<:AbstractFloat}
    Rₚ .= Rᵢ
    atom = sample(1:n_atoms(system), system.dynamics.moveable_atoms)
    apply_random_perturbation!(system, Rₚ, atom)
end

"""
    propose_centroid_move!(system::RingPolymerSystem{MonteCarlo}, Rᵢ::Array{T, 3}, Rₚ::Array{T, 3}) where {T<:AbstractFloat}
    
Propose a move for the ring polymer centroid for one atom.
"""
function propose_centroid_move!(system::RingPolymerSystem{MonteCarlo}, Rᵢ::Array{T, 3}, Rₚ::Array{T, 3}) where {T<:AbstractFloat}
    Rₚ .= Rᵢ
    transform_to_normal_modes!(system.ring_polymer, Rₚ, n_DoF(system))
    atom = sample(1:n_atoms(system), system.dynamics.moveable_atoms)
    @views apply_random_perturbation!(system, Rₚ[:,:,1], atom) # Perturb the centroid mode
    @views if atom ∉ system.ring_polymer.quantum_atoms # Ensure replicas are indentical if not quantum
        for i=2:n_beads(system)
            Rₚ[:, atom, i] .= Rₚ[:, atom, 1]
        end
    end
    transform_from_normal_modes!(system.ring_polymer, Rₚ, n_DoF(system))
end

"""
    apply_random_perturbation!(system::AbstractSystem{MonteCarlo}, R::AbstractMatrix, atom::Integer)
    
Randomly perturb the xyz coordinates of a single atom.
"""
function apply_random_perturbation!(system::AbstractSystem{MonteCarlo}, R::AbstractMatrix, atom::Integer)
    R[:, atom] .+= (rand(n_DoF(system)) .- 0.5) .* system.dynamics.Δ[system.atomic_parameters.atom_types[atom]] / sqrt(n_DoF(system))
end

"""
    propose_normal_mode_move!(system::RingPolymerSystem{MonteCarlo}, Rᵢ::Array{T, 3}, Rₚ::Array{T, 3}) where {T}
    
Propose a move for a single normal mode for a single atom.
"""
function propose_normal_mode_move!(system::RingPolymerSystem{MonteCarlo}, Rᵢ::Array{T, 3}, Rₚ::Array{T, 3}) where {T}
    Rₚ .= Rᵢ
    transform_to_normal_modes!(system.ring_polymer, Rₚ, n_DoF(system))
    atom = sample(system.ring_polymer.quantum_atoms)
    mode = sample(2:n_beads(system)-1)
    @views sample_mode!(system.ring_polymer, masses(system)[atom], mode, Rₚ[:, atom, mode])
    transform_from_normal_modes!(system.ring_polymer, Rₚ, n_DoF(system))
end

function sample_mode!(ring_polymer::Systems.RingPolymerParameters, mass::AbstractFloat, mode::Integer, R::AbstractArray)
    σ = sqrt(ring_polymer.ω_n / mass) / (ring_polymer.normal_mode_springs[mode] * mass)
    rand!(Normal(0, σ), R)
end

"""
    assess_proposal!(system::AbstractSystem{MonteCarlo}, Rᵢ, Rₚ, output, i)
    
Update the energy, check for acceptance, and update the output. 
"""
function assess_proposal!(system::AbstractSystem{MonteCarlo}, Rᵢ, Rₚ, output, i)
    system.dynamics.Eₚ = evaluate_configurational_energy(system, Rₚ)
    if acceptance_probability(system) > rand()
        system.dynamics.Eᵢ = system.dynamics.Eₚ
        Rᵢ .= Rₚ
        output.acceptance += 1
    end
    write_output!(output, Rᵢ, system.dynamics.Eᵢ, i)
end

"""
    apply_cell_boundaries!(system::System{MonteCarlo}, R::Matrix)
    
Ensure the atom remains inside the simulation cell.

Works only for cubic cells.
"""
function apply_cell_boundaries!(system::System{MonteCarlo}, R::Matrix)
    apply_cell_boundaries!(system.atomic_parameters.cell, R, n_atoms(system), n_DoF(system))
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
function apply_cell_boundaries!(system::RingPolymerSystem{MonteCarlo}, R::Array{T, 3}) where {T}
    transform_to_normal_modes!(system.ring_polymer, R, n_DoF(system))
    R[:,system.ring_polymer.quantum_atoms,1] ./= sqrt(n_beads(system))
    @views apply_cell_boundaries!(system.atomic_parameters.cell, R[:,:,1], n_atoms(system), n_DoF(system))
    R[:,system.ring_polymer.quantum_atoms,1] .*= sqrt(n_beads(system))
    transform_from_normal_modes!(system.ring_polymer, R, n_DoF(system))
end

"""
    apply_cell_boundaries!(cell::PeriodicCell, R::Matrix, n_atoms::Integer, n_DoF::Integer)
    
Apply simple periodic boundaries 

This works for cubic cells only.
"""
function apply_cell_boundaries!(cell::PeriodicCell, R::AbstractMatrix, n_atoms::Integer, n_DoF::Integer)
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
apply_cell_boundaries!(::InfiniteCell, ::AbstractMatrix, ::Integer, ::Integer) = nothing

"""
    acceptance_probability(system::AbstractSystem{MonteCarlo})

Return the Metropolis-Hastings acceptance probability.
"""
function acceptance_probability(system::AbstractSystem{MonteCarlo})
    acceptance_probability(system.dynamics.Eₚ, system.dynamics.Eᵢ, system.temperature)
end
function acceptance_probability(system::RingPolymerSystem{MonteCarlo})
    acceptance_probability(system.dynamics.Eₚ, system.dynamics.Eᵢ, system.temperature*n_beads(system))
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