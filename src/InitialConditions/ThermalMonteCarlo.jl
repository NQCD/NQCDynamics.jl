
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
  Estimators,
  DynamicsMethods,
  ndofs,
  nbeads,
  natoms,
  masses
using NQCCalculators

"""
    run_advancedmh_sampling(sim, r, steps, σ; movement_ratio=nothing, movement_ratio_internal=nothing, kwargs...)

Sample the configuration space for the simulation `sim` starting from `r`.

Total number of steps is given by `steps` and `σ` is the dictionary of
step sizes for each species.

`movement_ratio` denotes the fraction of system moved each step.
`internal_ratio` works as for `movement_ratio` but for the internal modes of the ring polymer.
For `movement_ratio = 0`, every degree of freedom is moved at each step, if `movement_ratio = 1`, then nothing will happen. 

If neither arguments are defined, default behaviour is to move one atom (and one ring polymer normal mode) per step on average. 

Further kwargs are passed to `AdvancedMH.sample` to allow for [extra functionality](https://turinglang.org/AbstractMCMC.jl/dev/api/#Common-keyword-arguments).
"""
function run_advancedmh_sampling(
  sim::AbstractSimulation,
  r,
  steps::Real,
  σ::Dict{Symbol,<:Real};
  movement_ratio=nothing,
  stop_ratio=nothing,
  movement_ratio_internal=nothing,
  stop_ratio_internal=nothing,
  move_ratio=nothing, # deprecated
  internal_ratio=nothing, # deprecated
  kwargs...
)

  # Give a warning if move_ratio or internal_ratio are used
  if move_ratio !== nothing || internal_ratio !== nothing
    @warn "move_ratio and internal_ratio kwargs are deprecated and may be removed in future. More information: https://nqcd.github.io/NQCDynamics.jl/stable/initialconditions/metropolishastings/"
    stop_ratio=move_ratio === nothing ? nothing : move_ratio
    stop_ratio_internal=internal_ratio === nothing ? nothing : internal_ratio
  end

  # If no kwargs for system fraction to move are given, perturb one atom and one normal mode at a time
  if movement_ratio===nothing && stop_ratio===nothing
    @debug "No movement restriction for atoms provided, automatically setting to move one atom per step on average."
    movement_ratio=1/size(sim)[2]
  end

  if movement_ratio_internal===nothing && stop_ratio_internal===nothing
    @debug "No movement restriction for ring polymer normal modes provided, automatically setting to move one mode per step on average."
    movement_ratio_internal=length(size(sim))==3 ? 1/size(sim)[3] : 0.0
  end

  # Set atom movement ratio by using whichever keyword is defined
  stop_ratio = stop_ratio===nothing ? 1-movement_ratio : stop_ratio
  stop_ratio_internal = stop_ratio_internal===nothing ? 1-movement_ratio_internal : stop_ratio_internal

  density = get_density_function(sim)

  density_model = AdvancedMH.DensityModel(density)
  proposal = get_proposal(sim, σ, stop_ratio, stop_ratio_internal)

  sampler = AdvancedMH.MetropolisHastings(proposal)

  initial_config = reshape_input(sim, copy(r))

  chain = AdvancedMH.sample(density_model, sampler, convert(Int, steps); initial_params=initial_config, kwargs...)

  return reshape_output(sim, chain)
end

function get_density_function(sim::Simulation)
  temperature = get_temperature(sim)
  function density(r)
    NQCCalculators.update_cache!(sim.cache, reshape(r, size(sim)))
    -DynamicsUtils.classical_potential_energy(sim, reshape(r, size(sim))) / temperature
  end
end

function get_deriv_function(sim::Simulation)
  temperature = get_temperature(sim)

  return function density_deriv(r)
    NQCCalculators.update_cache!(sim.cache, reshape(r, size(sim)))
    lπ = -DynamicsUtils.classical_potential_energy(sim, reshape(r, size(sim))) / temperature
    NQCCalculators.evaluate_derivative!(sim.cache, reshape(r, size(sim)))
    lπ_grad = -sim.cache.derivative / temperature
    return (lπ, lπ_grad[:])
  end

end

function get_density_function(sim::RingPolymerSimulation)
  temperature = get_ring_polymer_temperature(sim)

  return function density(r)
    #r_local = reshape(r, size(sim))
    NQCCalculators.update_cache!(sim.cache, reshape(r, size(sim)))
    RingPolymerArrays.transform_from_normal_modes!(reshape(r, size(sim)), sim.beads.transformation)
    potential = DynamicsUtils.classical_potential_energy(sim, reshape(r, size(sim)))
    spring = DynamicsUtils.classical_spring_energy(sim, reshape(r, size(sim)))
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
  stop_ratio::T
  sim::S
end

struct RingPolymerProposal{P,T,S<:RingPolymerSimulation} <: AdvancedMH.Proposal{P}
  proposal::P
  stop_ratio::T
  stop_ratio_internal::T
  sim::S
end

const MolecularProposal{P,T,S} = Union{RingPolymerProposal{P,T,S},ClassicalProposal{P,T,S}}

function ClassicalProposal(sim::Simulation, σ, stop_ratio)
  proposals = Matrix{UnivariateDistribution}(undef, size(sim))
  for (i, symbol) in enumerate(sim.atoms.types) # Position proposals
    distribution = σ[symbol] == 0 ? Dirac(0) : Normal(0, σ[symbol])
    proposals[:, i] .= distribution
  end
  ClassicalProposal(proposals[:], stop_ratio, sim)
end

function RingPolymerProposal(sim::RingPolymerSimulation, σ, stop_ratio, stop_ratio_internal)
  ωₖ = RingPolymers.get_matsubara_frequencies(nbeads(sim), sim.beads.ω_n)
  proposals = Array{UnivariateDistribution}(undef, size(sim))
  for i = 1:nbeads(sim)
    if i == 1
      for (j, symbol) in enumerate(sim.atoms.types)
        distribution = σ[symbol] == 0 ? Dirac(0.0) : Normal(0, σ[symbol])
        proposals[:, j, i] .= distribution
      end
    else
      for j in range(sim.atoms)
        if j in sim.beads.quantum_atoms
          Δ = sqrt(get_ring_polymer_temperature(sim) / (sim.atoms.masses[j] * ωₖ[i]^2))
          distribution = Normal(0, Δ)
        else
          distribution = Dirac(0)
        end
        proposals[:, j, i] .= distribution
      end
    end
  end
  RingPolymerProposal(proposals[:], stop_ratio, stop_ratio_internal, sim)
end

function Base.rand(rng::Random.AbstractRNG, p::ClassicalProposal)
  result = map(x -> rand(rng, x), p.proposal)
  result[Random.randsubseq(findall(x -> x!=0.0, result), p.stop_ratio)] .=0
  return result
end

function Base.rand(rng::Random.AbstractRNG, p::RingPolymerProposal)
  result = map(x -> rand(rng, x), p.proposal)
  reshaped_result = reshape(result, size(p.sim))

  # Zero some of the centroid moves
  sequence = Random.randsubseq(findall(x -> x != 0, CartesianIndices(size(p.sim)[1:2])), p.stop_ratio)
  reshaped_result[sequence, 1] .= 0

  # Zero some of the internal mode moves
  for i = 2:nbeads(p.sim)
    sequence = Random.randsubseq(findall(x -> x != 0, CartesianIndices((ndofs(p.sim), natoms(p.sim)))), p.stop_ratio_internal)
    reshaped_result[sequence, i] .= 0
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
  # proposal = AdvancedHMC.NUTS{AdvancedHMC.MultinomialTS, AdvancedHMC.GeneralisedNoUTurn}(integrator)
  kernel = AdvancedHMC.HMCKernel(AdvancedHMC.Trajectory{AdvancedHMC.MultinomialTS}(integrator, AdvancedHMC.GeneralisedNoUTurn()))
  adaptor = AdvancedHMC.StanHMCAdaptor(
    AdvancedHMC.MassMatrixAdaptor(metric),
    AdvancedHMC.StepSizeAdaptor(target_acceptance, integrator)
  )

  # Run the sampler to draw samples from the specified Gaussian, where
  #   - `samples` will store the samples
  #   - `stats` will store diagnostic statistics for each sample
  samples, stats = AdvancedHMC.sample(hamiltonian, kernel, r0, convert(Int, n_samples), adaptor; kwargs...)

  return reshape.(samples, Ref(size(sim))), stats
end

end # module
