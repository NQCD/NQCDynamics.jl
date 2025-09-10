
using NQCDistributions: DynamicalDistribution, ProductDistribution
using NQCCalculators: update_cache!

abstract type AbstractSelection end

"""
Select the initial conditions from the distribution in order. 
"""
struct OrderedSelection{D,I} <: AbstractSelection
    "Distribution that is sampled."
    distribution::D
    indices::I
end

"""
Obtain initial conditions by randomly sampling the distribution.
"""
struct RandomSelection{D} <: AbstractSelection
    "Distribution that is sampled."
    distribution::D
end

function Selection(distribution, selection::AbstractVector, trajectories)
    if trajectories > 1
        @info "Sampling the provided distribution in the range $selection." 
    end
    OrderedSelection(distribution, selection)
end

function Selection(distribution, ::Any, trajectories)
    if trajectories > 1
        @info "Sampling randomly from provided distribution." 
    end
    RandomSelection(distribution)
end

function (select::OrderedSelection)(prob, i, repeat)
    j = select.indices[i]
    u0 = sample_distribution(prob.p, select.distribution, j)
    DynamicsMethods.create_problem(u0, prob.tspan, prob.p)
end

"""
    sample_distribution(sim::AbstractSimulation, distribution::DynamicalDistribution, i)

Generate initial conditions for a dynamics simulation as a selection, `i`, from the provided distribution of nuclear degrees of freedom.
Simulation cache is updated to according to the sampled positions.\n
A container of dynamics variables are created from the intial condtions and returned in the format of a `ComponentVector`.
"""
function sample_distribution(sim::AbstractSimulation, distribution::DynamicalDistribution, i)
    u = distribution[i]
    update_cache!(sim.cache, copy(u.r))
    DynamicsMethods.DynamicsVariables(sim, u.v, u.r)
end

"""
    sample_distribution(sim::AbstractSimulation, distribution::ProductDistribution, i)

Generate initial conditions for a dynamics simulation as a selection, `i`, from the provided product distribution of nuclear degrees of freedom and the electronic state.
Simulation cache is updated to according to the sampled positions.\n
A container of dynamics variables are created from the intial condtions and returned in the format of a `ComponentVector`.
"""
function sample_distribution(sim::AbstractSimulation, distribution::ProductDistribution, i)
    u = distribution.nuclear[i]
    update_cache!(sim.cache, copy(u.r))
    DynamicsMethods.DynamicsVariables(sim, copy(u.v), copy(u.r), distribution.electronic)
end

function (select::RandomSelection)(prob, i, repeat)
    u0 = sample_distribution(prob.p, select.distribution)
    DynamicsMethods.create_problem(u0, prob.tspan, prob.p)
end

"""
    sample_distribution(sim::AbstractSimulation, distribution::DynamicalDistribution)

Generate initial conditions for a dynamics simulation by randomly sampling from the provided distribution of nuclear degrees of freedom.
Simulation cache is updated to according to the sampled positions.\n
A container of dynamics variables are created from the intial condtions and returned in the format of a `ComponentVector`.
"""
function sample_distribution(sim::AbstractSimulation, distribution::DynamicalDistribution)
    u = rand(distribution)
    update_cache!(sim.cache, copy(u.r))
    DynamicsMethods.DynamicsVariables(sim, u.v, u.r)
end

"""
    sample_distribution(sim::AbstractSimulation, distribution::ProductDistribution)

Generate initial conditions for a dynamics simulation by randomly sampling from the provided product distribution of nuclear degrees of freedom and the electronic state.
Simulation cache is updated to according to the sampled positions.\n
A container of dynamics variables are created from the intial condtions and returned in the format of a `ComponentVector`.
"""
function sample_distribution(sim::AbstractSimulation, distribution::ProductDistribution)
    u = rand(distribution.nuclear)
    update_cache!(sim.cache, copy(u.r))
    DynamicsMethods.DynamicsVariables(sim, copy(u.v), copy(u.r), distribution.electronic)
end

sample_distribution(::AbstractSimulation, distribution) = copy(distribution)

