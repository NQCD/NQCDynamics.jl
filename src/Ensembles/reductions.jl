using StatsBase: mean
using NonadiabaticMolecularDynamics: TimeCorrelationFunctions
using .TimeCorrelationFunctions: TimeCorrelationFunction

abstract type AbstractReduction end

"""
Sum the outputs from each trajectory.
"""
struct SumReduction <: AbstractReduction end
(reduction::SumReduction)(u, batch, I) = u + sum(batch), false

"""
Average the outputs over all trajectories.
"""
mutable struct MeanReduction <: AbstractReduction
    counter::Int
end

function (reduction::MeanReduction)(u, batch, I)
    failed_trajectories = length.(batch) .!= length(u)
    deleteat!(batch, failed_trajectories)

    previous = u .* reduction.counter
    reduction.counter += length(batch)
    current = previous .+ sum(batch)

    (current ./ reduction.counter, false)
end

function select_reduction(reduction::Symbol)
    if reduction === :mean
        return MeanReduction(0)
    elseif reduction === :sum
        return SumReduction()
    elseif reduction === :append
        return (u,data,I)->(append!(u,data),false)
    else
        throw(ArgumentError("`reduction` $reduction not recognised."))
    end
end

function get_u_init(::Union{MeanReduction, SumReduction}, saveat, stripped_kwargs, tspan, u0, output::TimeCorrelationFunction)
    if saveat != []
        if saveat isa Number
            savepoints = tspan[1]:saveat:tspan[2]
        else
            savepoints = saveat
        end
    elseif haskey(stripped_kwargs, :dt)
        dt = stripped_kwargs[:dt]
        savepoints = tspan[1]:dt:tspan[2]
    else
        throw(ArgumentError("Make sure to pass either `saveat` or `dt` as keyword arguments
            to compute average/sum over many trajectories."))
    end
    u_init = [TimeCorrelationFunctions.correlation_template(output) for _ ∈ savepoints]
    return u_init
end

function get_u_init(::Union{MeanReduction, SumReduction}, saveat, stripped_kwargs, tspan, u0, output::AbstractOutput)
    return output_template(output, u0)
end

get_u_init(reduction, saveat, stripped_kwargs, tspan, u0, output) = []
