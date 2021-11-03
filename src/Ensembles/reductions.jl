using StatsBase: mean

abstract type AbstractReduction end

"""
Sum the outputs from each trajectory.
"""
struct SumReduction{T} <: AbstractReduction
    "Template for the reduction variable."
    u_init::T
end
SumReduction() = SumReduction(0.0)

(reduction::SumReduction)(u, batch, I) = u + sum(batch), false

"""
Average the outputs over all trajectories.
"""
mutable struct MeanReduction{T} <: AbstractReduction
    u_init::T
    counter::Int
end
MeanReduction(u_init) = MeanReduction(u_init, 0)
MeanReduction() = MeanReduction(0.0)

function (reduction::MeanReduction)(u, batch, I)
    failed_trajectories = length.(batch) .!= length(reduction.u_init)
    deleteat!(batch, failed_trajectories)

    previous = u .* reduction.counter
    reduction.counter += length(batch)
    current = previous .+ sum(batch)

    (current ./ reduction.counter, false)
end
