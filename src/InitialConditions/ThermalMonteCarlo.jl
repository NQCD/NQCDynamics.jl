
module ThermalMonteCarlo

using Distributions: UnivariateDistribution, Dirac, Normal
using ComponentArrays: ComponentVector
using AdvancedMH: AdvancedMH
using AdvancedHMC: AdvancedHMC
using Random: Random
using UnPack: @unpack
using RingPolymerArrays: RingPolymerArrays

using NQCDynamics:
    NQCDynamics,
    AbstractSimulation,
    Simulation,
    RingPolymerSimulation,
    get_temperature,
    get_ring_polymer_temperature,
    DynamicsUtils,
    RingPolymers,
    Calculators,
    Estimators,
    DynamicsMethods,
    ndofs,
    nbeads,
    natoms,
    masses

"""
    run_advancedhmc_sampling(sim, r, steps, σ; move_ratio=0.0, internal_ratio=0.0)

Sample the configuration space for the simulation `sim` starting from `r`.

Total number of steps is given by `steps` and `σ` is the dictionary of
step sizes for each species.

`move_ratio` defaults to `0.0` and denotes the fraction of system moved each step.
If `move_ratio = 0`, every degree of freedom is moved at each step.
If `move_ratio = 1`, then nothing will happen. Experiment with this parameter to achieve
optimal sampling.

`internal_ratio` works as for `move_ratio` but for the internal modes of the ring polymer.
"""
function run_advancedmh_sampling(
    sim::AbstractSimulation,
    r,
    steps::Real,
    σ::Dict{Symbol,<:Real};
    move_ratio=0.0,
    internal_ratio=0.0
    )

    density = get_density_function(sim)

    density_model = AdvancedMH.DensityModel(density)
    proposal = get_proposal(sim, σ, move_ratio, internal_ratio)

    sampler = AdvancedMH.MetropolisHastings(proposal)

    initial_config = reshape_input(sim, copy(r))

    chain = AdvancedMH.sample(density_model, sampler, convert(Int, steps);
                              init_params=initial_config)

    return reshape_output(sim, chain)
end

function get_density_function(sim::Simulation)
    temperature = get_temperature(sim)
    function density(r)
        -DynamicsUtils.classical_potential_energy(sim, reshape(r, size(sim))) / temperature
    end
end

function get_deriv_function(sim::Simulation)
    temperature = get_temperature(sim)
    function density_deriv(r)
        lπ = -DynamicsUtils.classical_potential_energy(sim, reshape(r, size(sim))) / temperature
        Calculators.evaluate_derivative!(sim.calculator, reshape(r, size(sim)))
        lπ_grad = -sim.calculator.derivative / temperature
        return (lπ, lπ_grad[:])
    end
end

function get_density_function(sim::RingPolymerSimulation)
    temperature = get_ring_polymer_temperature(sim)
    function density(r)
        r_local = reshape(copy(r), size(sim))
        RingPolymerArrays.transform_from_normal_modes!(r_local, sim.beads.transformation)
        potential = DynamicsUtils.classical_potential_energy(sim, r_local)
        spring = DynamicsUtils.classical_spring_energy(sim, r_local)
        total = potential + spring
        return -total / temperature
    end
end

get_proposal(sim::Simulation, σ, move_ratio, internal_ratio) = ClassicalProposal(sim, σ, move_ratio)
get_proposal(sim::RingPolymerSimulation, σ, move_ratio, internal_ratio) = RingPolymerProposal(sim, σ, move_ratio, internal_ratio)

reshape_input(::Simulation, r) = r[:]
function reshape_input(sim::RingPolymerSimulation, r)
    RingPolymerArrays.transform_to_normal_modes!(r, sim.beads.transformation)
    return r[:]
end

function reshape_output(sim::Simulation, chain)
    [reshape(config.params, size(sim)) for config in chain]
end
function reshape_output(sim::RingPolymerSimulation, chain)
    rs = [reshape(copy(config.params), size(sim)) for config in chain]
    for r in rs
        RingPolymerArrays.transform_from_normal_modes!(r, sim.beads.transformation)
    end
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
    internal_ratio::T
    sim::S
end

const MolecularProposal{P,T,S} = Union{RingPolymerProposal{P,T,S}, ClassicalProposal{P,T,S}}

function ClassicalProposal(sim::Simulation, σ, move_ratio)
    proposals = Matrix{UnivariateDistribution}(undef, size(sim))
    for (i, symbol) in enumerate(sim.atoms.types) # Position proposals
        distribution = σ[symbol] == 0 ? Dirac(0) : Normal(0, σ[symbol])
        proposals[:,i] .= distribution
    end
    ClassicalProposal(proposals[:], move_ratio, sim)
end

function RingPolymerProposal(sim::RingPolymerSimulation, σ, move_ratio, internal_ratio)
    ωₖ = RingPolymers.get_matsubara_frequencies(nbeads(sim), sim.beads.ω_n)
    proposals = Array{UnivariateDistribution}(undef, size(sim))
    for i=1:nbeads(sim)
        if i == 1
            for (j, symbol) in enumerate(sim.atoms.types)
                distribution = σ[symbol] == 0 ? Dirac(0.0) : Normal(0, σ[symbol])
                proposals[:,j,i] .= distribution
            end
        else
            for j in range(sim.atoms)
                if j in sim.beads.quantum_atoms
                    Δ = sqrt(get_ring_polymer_temperature(sim) / (sim.atoms.masses[j] * ωₖ[i]^2))
                    distribution = Normal(0, Δ)
                else
                    distribution = Dirac(0)
                end
                proposals[:,j,i] .= distribution
            end
        end
    end
    RingPolymerProposal(proposals[:], move_ratio, internal_ratio, sim)
end

function Base.rand(rng::Random.AbstractRNG, p::ClassicalProposal)
    result = map(x -> rand(rng, x), p.proposal)
    result[Random.randsubseq(eachindex(result), p.move_ratio)] .= 0
    return result
end

function Base.rand(rng::Random.AbstractRNG, p::RingPolymerProposal)
    result = map(x -> rand(rng, x), p.proposal)
    reshaped_result = reshape(result, size(p.sim))

    # Zero some of the centroid moves
    sequence = Random.randsubseq(CartesianIndices(size(p.sim)[1:2]), p.move_ratio)
    reshaped_result[sequence,1] .= 0

    # Zero some of the internal mode moves
    for i=2:nbeads(p.sim)
        sequence = Random.randsubseq(CartesianIndices((ndofs(p.sim), natoms(p.sim))), p.internal_ratio)
        reshaped_result[sequence,i] .= 0
    end

    return result
end

function AdvancedMH.propose(rng::Random.AbstractRNG, proposal::MolecularProposal,
                            ::AdvancedMH.DensityModel)
    return rand(rng, proposal)
end

function AdvancedMH.propose(rng::Random.AbstractRNG, proposal::MolecularProposal,
                            ::AdvancedMH.DensityModel, t)
    return t + rand(rng, proposal)
end

AdvancedMH.logratio_proposal_density(::MolecularProposal, state, candidiate) = 0

"""
    run_advancedhmc_sampling(sim::Simulation, r, n_samples;
        target_acceptance=0.5, kwargs...)

Perform Hamiltonian Monte Carlo sampling for the simulation `sim` using `AdvancedHMC.jl`.
"""
function run_advancedhmc_sampling(sim::Simulation, r, n_samples; target_acceptance=0.5, kwargs...)

    density = get_density_function(sim)
    deriv = get_deriv_function(sim)

    r0 = reshape_input(sim, copy(r))

    # Define a Hamiltonian system
    metric = AdvancedHMC.DiagEuclideanMetric(length(r0))
    hamiltonian = AdvancedHMC.Hamiltonian(metric, density, deriv)

    # Define a leapfrog solver, with initial step size chosen heuristically
    initial_ϵ = AdvancedHMC.find_good_stepsize(hamiltonian, r0)
    integrator = AdvancedHMC.Leapfrog(initial_ϵ)

    # Define an HMC sampler, with the following components
    #   - multinomial sampling scheme,
    #   - generalised No-U-Turn criteria, and
    #   - windowed adaption for step-size and diagonal mass matrix
    proposal = AdvancedHMC.NUTS{AdvancedHMC.MultinomialTS, AdvancedHMC.GeneralisedNoUTurn}(integrator)
    adaptor = AdvancedHMC.StanHMCAdaptor(
        AdvancedHMC.MassMatrixAdaptor(metric),
        AdvancedHMC.StepSizeAdaptor(target_acceptance, integrator)
        )

    # Run the sampler to draw samples from the specified Gaussian, where
    #   - `samples` will store the samples
    #   - `stats` will store diagnostic statistics for each sample
    samples, stats = AdvancedHMC.sample(hamiltonian, proposal, r0, convert(Int, n_samples), adaptor; kwargs...)

    return reshape.(samples, Ref(size(sim))), stats
end

end # module
