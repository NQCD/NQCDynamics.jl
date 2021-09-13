import AdvancedMH
using AdvancedMH: DensityModel, MetropolisHastings, RandomWalkProposal, sample
using Distributions: Dirac
using Random
using ComponentArrays: ComponentVector

function sample_configurations(
    sim::AbstractSimulation,
    u0,
    steps::Real,
    σ::Dict{Symbol,<:Real};
    move_ratio=0.9
    )

    shape = size(get_positions(u0))
    density = get_density_function(sim, shape)

    density_model = DensityModel(density)
    proposal = get_proposal(sim, σ, move_ratio)

    sampler = AdvancedMH.MetropolisHastings(proposal)

    initial_config = reshape_input(sim, copy(u0))

    chain = sample(density_model, sampler, convert(Int, steps); init_params=initial_config)

    return reshape_output(sim, chain, shape)
end

function get_density_function(sim::Simulation, shape)
    temperature = get_temperature(sim)
    function density(u)
        v = @view u[1:prod(shape)]
        r = @view u[prod(shape)+1:end]
        -evaluate_hamiltonian(sim, reshape(v, shape), reshape(r, shape)) / temperature
    end
end

function get_density_function(sim::RingPolymerSimulation, shape)
    temperature = get_temperature(sim)
    function density(u)
        v = @view u[1:prod(shape)]
        r = @view u[prod(shape)+1:end]
        v = RingPolymerArray(reshape(copy(v), shape); normal=true)
        r = RingPolymerArray(reshape(copy(r), shape); normal=true)
        -evaluate_hamiltonian(sim, v, r) / temperature
    end
end

get_proposal(sim::Simulation, σ, move_ratio) =  ClassicalProposal(sim, σ, 1-move_ratio)
get_proposal(sim::RingPolymerSimulation, σ, move_ratio) = RingPolymerProposal(sim, σ, 1-move_ratio)

reshape_input(::Simulation, u0) = u0[:]
function reshape_input(sim::RingPolymerSimulation, u0)
    r = get_positions(u0)
    v = get_velocities(u0)
    transform_to_normal_modes!(sim.beads, r)
    transform_to_normal_modes!(sim.beads, v)
    return vcat(v[:], r[:])
end

function reshape_output(::Simulation, chain, shape)
    s = prod(shape)
    [ComponentVector(v=reshape(config.params[1:s], shape), r=reshape(config.params[s+1:end], shape)) for config in chain]
end
function reshape_output(sim::RingPolymerSimulation, chain, shape)
    s = prod(shape)
    rs = [ComponentVector(v=reshape(config.params[1:s], shape), r=reshape(config.params[s+1:end], shape)) for config in chain]
    transform_from_normal_modes!.(sim.beads, rs)
    return rs
end

struct ClassicalProposal{P,T,S<:Simulation} <: AdvancedMH.Proposal{P}
    proposal::P
    move_ratio::T
    sim::S
end

struct RingPolymerProposal{P,T,S<:RingPolymerSimulation} <: AdvancedMH.Proposal{P}
    proposal::P
    move_ratio::T
    sim::S
end

const MolecularProposal{P,T,S} = Union{RingPolymerProposal{P,T,S}, ClassicalProposal{P,T,S}}

function ClassicalProposal(sim, σ, move_ratio)
    proposals = Array{Distributions.UnivariateDistribution,3}(undef, sim.DoFs, length(sim.atoms), 2)
    for i in range(sim.atoms) # Velocity proposals
        distribution = Normal(0, sqrt(get_temperature(sim) / sim.atoms.masses[i]))
        proposals[:,i,1] .= distribution
    end
    for (i, symbol) in enumerate(sim.atoms.types) # Position proposals
        distribution = σ[symbol] == 0 ? Dirac(0) : Normal(0, σ[symbol])
        proposals[:,i,2] .= distribution
    end
    ClassicalProposal(proposals[:], move_ratio, sim)
end

function RingPolymerProposal(sim, σ, move_ratio)
    ωₖ = NonadiabaticMolecularDynamics.get_matsubara_frequencies(length(sim.beads), sim.beads.ω_n)
    proposals = Array{Distributions.UnivariateDistribution,4}(undef, sim.DoFs, length(sim.atoms), length(sim.beads), 2)
    for i=1:length(sim.beads)
        if i == 1
            for j in range(sim.atoms)
                distribution = Normal(0, sqrt(get_temperature(sim) / sim.atoms.masses[j]))
                proposals[:,j,i,1] .= distribution
            end
        else
            for j in range(sim.atoms)
                if j in sim.beads.quantum_atoms
                    Δ = sqrt(get_temperature(sim) / sim.atoms.masses[j])
                    distribution = Normal(0, Δ)
                else
                    distribution = Dirac(0.0)
                end
                proposals[:,j,i,1] .= distribution
            end
        end
    end
    for i=1:length(sim.beads)
        if i == 1
            for (j, symbol) in enumerate(sim.atoms.types)
                distribution = σ[symbol] == 0 ? Dirac(0.0) : Normal(0, σ[symbol])
                proposals[:,j,i,2] .= distribution
            end
        else
            for j in range(sim.atoms)
                if j in sim.beads.quantum_atoms
                    Δ = sqrt(get_temperature(sim) / (sim.atoms.masses[j] * ωₖ[i]^2))
                    distribution = Normal(0, Δ)
                else
                    distribution = Dirac(0)
                end
                proposals[:,j,i,2] .= distribution
            end
        end
    end
    RingPolymerProposal(proposals[:], move_ratio, sim)
end

function Base.rand(rng::Random.AbstractRNG, p::MolecularProposal)
    result = map(x -> rand(rng, x), p.proposal)
    result[randsubseq(eachindex(result), p.move_ratio)] .= 0
    return result
end

function AdvancedMH.propose(rng::Random.AbstractRNG, proposal::MolecularProposal, ::DensityModel)
    return rand(rng, proposal)
end

function AdvancedMH.propose(rng::Random.AbstractRNG, proposal::MolecularProposal, ::DensityModel, t)
    return t + rand(rng, proposal)
end

AdvancedMH.logratio_proposal_density(::MolecularProposal, state, candidiate) = 0
