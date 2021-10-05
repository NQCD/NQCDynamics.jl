using StatsBase: mean

abstract type AbstractReduction end

"""
Sum the outputs from each trajectory.
"""
struct SumReduction{T} <: AbstractReduction
    "Template for the reduction variable."
    u_init::T
    SumReduction(u_init) = new{typeof(u_init)}(u_init)
end
SumReduction() = SumReduction(0.0)

(reduction::SumReduction)(u, batch, I) = u + sum(batch), false

"""
Average the outputs over all trajectories.
"""
struct MeanReduction <: AbstractReduction end
(reduction::MeanReduction)(u, batch, I) = (mean(batch), false)
