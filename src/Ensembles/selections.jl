using DiffEqBase: CallbackSet

using NQCDynamics: InitialConditions, DynamicsOutputs

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

function Selection(distribution, selection::AbstractVector)
    @info "Sampling the provided distribution in the range $selection." 
    OrderedSelection(distribution, selection)
end

function Selection(distribution, ::Any)
    @info "Sampling randomly from provided distribution." 
    RandomSelection(distribution)
end

function (select::OrderedSelection)(prob, i, repeat)
    j = select.indices[i]
    u0 = sample_distribution(prob.p, select.distribution, j)
    DynamicsMethods.create_problem(u0, prob.tspan, prob.p)
end

function (select::RandomSelection)(prob, i, repeat)
    j = rand(1:lastindex(select.distribution))
    u0 = sample_distribution(prob.p, select.distribution, j)
    DynamicsMethods.create_problem(u0, prob.tspan, prob.p)
end
