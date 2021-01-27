
export PositionDistribution
export PhasespaceDistribution

using Random: AbstractRNG, rand!
using Distributions

struct PhasespaceVariate <: VariateForm end

struct PhasespaceDistribution{R,P,S} <: Sampleable{PhasespaceVariate,Continuous}
    positions::R
    momenta::P
    size::NTuple{S,Int}
end

Base.eltype(s::PhasespaceDistribution{<:Sampleable,P}) where {P} = eltype(s.positions)
Base.eltype(s::PhasespaceDistribution{<:AbstractArray,P} where {P}) = eltype(s.positions[1])
Base.size(s::PhasespaceDistribution) = s.size

function Distributions.rand(rng::AbstractRNG, s::Sampleable{PhasespaceVariate})
    Distributions._rand!(rng, s, [Array{eltype(s)}(undef, size(s)) for i=1:2])
end

function Distributions._rand!(rng::AbstractRNG, s::PhasespaceDistribution, x::Vector{<:Array})
    i = rand(rng, 1:length(s.positions))
    x[1] .= select_item(s.positions, i, s.size)
    x[2] .= select_item(s.momenta, i, s.size)
    x
end

select_item(x::Vector, i::Integer, ::NTuple) = x[i]
select_item(x::Sampleable{Univariate}, ::Integer, size::NTuple) = rand(x, size)
select_item(x::Real, ::Integer, ::NTuple) = x
