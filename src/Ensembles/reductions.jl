using StatsBase: mean

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
