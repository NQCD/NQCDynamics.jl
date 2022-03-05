using StatsBase: mean
using NQCDynamics: TimeCorrelationFunctions, DynamicsOutputs
using .TimeCorrelationFunctions: TimeCorrelationFunction

"""
Sum the outputs from each trajectory.
"""
struct SumReduction end
function (reduction::SumReduction)(u, batch, I)
    sum_outputs!(u, batch)
    return (u, false)
end

"""
Average the outputs over all trajectories.
"""
mutable struct MeanReduction
    counter::Int
end

function (reduction::MeanReduction)(u, batch, I)
    u .*= reduction.counter
    reduction.counter += length(batch)
    sum_outputs!(u, batch)
    u ./= reduction.counter
    return (u, false)
end

function sum_outputs!(u, batch)
    for i=1:length(batch)
        if size(batch[i]) != size(u)
            @warn "DimensionMismatch encountered when reducing outputs. \
                Discarding trajectory. Check `u_init` matches the output size."
        else
            u .+= batch[i]
        end
    end
end

struct AppendReduction end
(::AppendReduction)(u,data,I) = (append!(u,data), false)

"""
    Reduction(reduction::Symbol)

Converts reduction keyword into function that performs reduction.

* `:discard` is used internally when using callbacks to save outputs instead.
"""
function Reduction(reduction::Symbol)
    if reduction === :mean
        return MeanReduction(0)
    elseif reduction === :sum
        return SumReduction()
    elseif reduction === :append
        return AppendReduction()
    else
        throw(ArgumentError("`reduction` $reduction not recognised."))
    end
end

function get_u_init(::Union{MeanReduction, SumReduction}, saveat, stripped_kwargs, tspan, u0, output::TimeCorrelationFunction)
    savepoints = get_savepoints(saveat, stripped_kwargs, tspan)
    u_init = [TimeCorrelationFunctions.correlation_template(output) for _ ∈ savepoints]
    return u_init
end

function get_u_init(::Union{MeanReduction, SumReduction}, saveat, stripped_kwargs, tspan, u0, output::DynamicsOutputs.EnsembleSaver)
    savepoints = get_savepoints(saveat, stripped_kwargs, tspan)
    u_init = DynamicsOutputs.output_template(output, savepoints, u0)
    return u_init
end

function get_u_init(::Union{MeanReduction, SumReduction}, saveat, stripped_kwargs, tspan, u0, output)
    return output_template(output, u0)
end

function get_u_init(::AppendReduction, saveat, stripped_kwargs, tspan, u0, output)
    return nothing
end

function get_savepoints(saveat, stripped_kwargs, tspan)
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
end
