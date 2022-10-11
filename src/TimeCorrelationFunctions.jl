
"""
    TimeCorrelationFunctions

This module defines extra types that can be used as Ensemble outputs when computing
time-correlation functions.
It hopes to provide a minimal interface that reduces code repetition when implementing
different correlation functions.
"""
module TimeCorrelationFunctions

using NQCDynamics: AbstractSimulation, Estimators
using NQCDynamics.NQCModels: nstates

using NQCDistributions: Diabatic, Adiabatic

export PopulationCorrelationFunction

"""
    TimeCorrelationFunction

Abstract type for defining time correlation functions
"""
abstract type TimeCorrelationFunction end

function (correlation::TimeCorrelationFunction)(sol, _)
    sim = correlation.sim
    template = correlation_template(correlation)
    out = [similar(template) for _ in 1:length(sol.u)]

    normalisation = evaluate_normalisation(sim, correlation)
    initial_value = evaluate_initial_value(sim, correlation, first(sol.u))

    for i in 1:length(sol.u)
        final_value = evaluate_final_value(sim, correlation, sol.u[i])
        correlate!(out[i], normalisation, initial_value, final_value, correlation)
    end

    return out
end

evaluate_normalisation(sim::AbstractSimulation, correlation::TimeCorrelationFunction) = 1
evaluate_initial_value(sim::AbstractSimulation, correlation::TimeCorrelationFunction, u) = 
    evaluate_final_value(sim, correlation, u)

function correlate!(out, normalisation, initial_value, final_value, ::TimeCorrelationFunction)
    @. out = normalisation * initial_value * final_value
end

function correlation_template(::TimeCorrelationFunction)
    error("Implement this for your correlation function type.")
end

"""
    PopulationCorrelationFunction{T,S<:AbstractSimulation} <: TimeCorrelationFunction

Output type for computing the population correlation function.
The `statetype` determines the population type (diabatic or adiabatic).
`sim` must also be provided to access the parameters to compute the population.
"""
struct PopulationCorrelationFunction{T,S<:AbstractSimulation} <: TimeCorrelationFunction
    sim::S
    statetype::T
end

function correlation_template(correlation::PopulationCorrelationFunction)
    zeros(nstates(correlation.sim), nstates(correlation.sim))
end

function evaluate_final_value(sim::AbstractSimulation, ::PopulationCorrelationFunction{Diabatic}, u)
    Estimators.diabatic_population(sim, u)
end

function evaluate_initial_value(sim::AbstractSimulation, ::PopulationCorrelationFunction{Diabatic}, u)
    Estimators.initial_diabatic_population(sim, u)
end

function evaluate_final_value(sim::AbstractSimulation, ::PopulationCorrelationFunction{Adiabatic}, u)
    Estimators.adiabatic_population(sim, u)
end

function evaluate_initial_value(sim::AbstractSimulation, ::PopulationCorrelationFunction{Adiabatic}, u)
    Estimators.initial_adiabatic_population(sim, u)
end

function correlate!(out, normalisation::Real, initial_value::Vector, final_value::Vector, ::PopulationCorrelationFunction)
    out .= normalisation .* initial_value * final_value'
end

end # module
