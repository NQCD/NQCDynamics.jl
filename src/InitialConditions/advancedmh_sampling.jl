import AdvancedMH
using AdvancedMH: DensityModel, MetropolisHastings, RandomWalkProposal, sample
using Distributions: MvNormal
using Random

function sample_configurations(
    sim::AbstractSimulation,
    R0::AbstractArray,
    steps::Real,
    σ::Dict{Symbol,T};
    kwargs...
    ) where {T}

    shape = size(R0)
    density = get_density_function(sim, shape)

    density_model = DensityModel(density)
    proposal = get_proposal(sim, σ; kwargs...)

    sampler = AdvancedMH.MetropolisHastings(proposal)

    initial_config = reshape_input(sim, copy(R0))

    chain = sample(density_model, sampler, convert(Int, steps); init_params=initial_config)

    return reshape_output(sim, chain, shape)
end

function get_density_function(sim::Simulation, shape)
    temperature = get_temperature(sim)
    density(R) = -evaluate_potential_energy(sim, reshape(R, shape)) / temperature
end

function get_density_function(sim::RingPolymerSimulation, shape)
    temperature = get_temperature(sim)
    function density(R)
        r = RingPolymerArray(reshape(copy(R), shape); normal=true)
        -evaluate_potential_energy(sim, r) / temperature
    end
end

function get_proposal(sim::Simulation, σ)
    proposals = Array{Distributions.UnivariateDistribution,2}(undef, sim.DoFs, length(sim.atoms))
    for (i, symbol) in enumerate(sim.atoms.types)
        distribution = σ[symbol] == 0 ? Dirac(0) : Normal(0, σ[symbol])
        for j=1:sim.DoFs
            proposals[j,i] = distribution
        end
    end
    return RandomWalkProposal(proposals[:])
end
function get_proposal(sim::RingPolymerSimulation, σ; move_ratio=0.9)
    return RingPolymerProposal(sim, σ, 1-move_ratio)
end

reshape_input(::Simulation, R0::AbstractMatrix) = R0[:]
function reshape_input(sim::RingPolymerSimulation, R0)
    transform_to_normal_modes!(sim.beads, R0)
    return R0[:]
end

reshape_output(::Simulation, chain, shape) = [reshape(config.params, shape) for config in chain]
function reshape_output(sim::RingPolymerSimulation, chain, shape)
    rs = [reshape(copy(config.params), shape) for config in chain]
    transform_from_normal_modes!.(sim.beads, rs)
    return rs
end

struct RingPolymerProposal{P,T,S<:RingPolymerSimulation} <: AdvancedMH.Proposal{P}
    proposal::P
    move_ratio::T
    sim::S
end

function RingPolymerProposal(sim, σ, move_ratio)
    ωₖ = NonadiabaticMolecularDynamics.get_matsubara_frequencies(length(sim.beads), sim.beads.ω_n)
    proposals = Array{Distributions.UnivariateDistribution,3}(undef, sim.DoFs, length(sim.atoms), length(sim.beads))
    for i=1:length(sim.beads)
        if i == 1
            for (j, symbol) in enumerate(sim.atoms.types)
                distribution = σ[symbol] == 0 ? Dirac(0.0) : Normal(0, σ[symbol])
                for k=1:sim.DoFs
                    proposals[k,j,i] = distribution
                end
            end
        else
            for j in sim.beads.quantum_atoms
                for k=1:sim.DoFs
                    Δ = sqrt(get_temperature(sim) / (sim.atoms.masses[j] * ωₖ[i]^2))
                    proposals[k,j,i] = Normal(0, Δ)
                end
            end
        end
    end
    RingPolymerProposal(proposals[:], move_ratio, sim)
end

function Base.rand(rng::Random.AbstractRNG, p::RingPolymerProposal{<:AbstractArray})
    result = map(x -> rand(rng, x), p.proposal)

    shaped = reshape(result, p.sim.DoFs, length(p.sim.atoms), length(p.sim.beads))
    for i in randsubseq(rng, 2:length(p.sim.beads), p.move_ratio)
        for j in p.sim.beads.quantum_atoms
            for k=1:p.sim.DoFs
                shaped[k,j,i] = 0
            end
        end
    end

    return result
end

function AdvancedMH.propose(rng::Random.AbstractRNG, proposal::RingPolymerProposal, ::DensityModel)
    return rand(rng, proposal)
end

function AdvancedMH.propose(rng::Random.AbstractRNG, proposal::RingPolymerProposal, ::DensityModel, t)
    return t + rand(rng, proposal)
end

function AdvancedMH.q(proposal::RingPolymerProposal, t, t_cond)
    return logpdf(proposal, t - t_cond)
end

AdvancedMH.logratio_proposal_density(::RingPolymerProposal, state, candidiate) = 0
