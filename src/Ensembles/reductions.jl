using StatsBase: mean
using NQCDynamics: TimeCorrelationFunctions
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

"""
    select_reduction(reduction::Symbol)

Converts reduction keyword into function that performs reduction.

* `:discard` is used internally when using callbacks to save outputs instead.
"""
function select_reduction(reduction::Symbol)
    if reduction === :mean
        return MeanReduction(0)
    elseif reduction === :sum
        return SumReduction()
    elseif reduction === :append
        return (u,data,I)->(append!(u,data),false)
    elseif reduction === :discard
        return (u,data,I)->((nothing,),false)
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

function get_u_init(::Union{MeanReduction, SumReduction}, saveat, stripped_kwargs, tspan, u0, output)
    return output_template(output, u0)
end

get_u_init(reduction, saveat, stripped_kwargs, tspan, u0, output) = Nothing[]
