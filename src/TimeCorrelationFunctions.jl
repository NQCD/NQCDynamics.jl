
module TimeCorrelationFunctions

using NonadiabaticMolecularDynamics: AbstractSimulation, Estimators
using NonadiabaticMolecularDynamics.NonadiabaticDistributions: Diabatic, Adiabatic
using NonadiabaticMolecularDynamics.NonadiabaticModels: nstates

export PopulationCorrelationFunction

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

    return (out, false)
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
