
export DynamicalDistribution

using Random: AbstractRNG, rand!
using Distributions

struct DynamicalVariate <: VariateForm end

struct DynamicalDistribution{V,R,S} <: Sampleable{DynamicalVariate,Continuous}
    velocities::V
    positions::R
    size::NTuple{S,Int}
end

Base.eltype(s::DynamicalDistribution{<:Sampleable,R}) where {R} = eltype(s.velocities)
Base.eltype(s::DynamicalDistribution{<:AbstractArray,R} where {R}) = eltype(s.velocities[1])
Base.size(s::DynamicalDistribution) = s.size

function Distributions.rand(rng::AbstractRNG, s::Sampleable{DynamicalVariate})
    Distributions._rand!(rng, s, [Array{eltype(s)}(undef, size(s)) for i=1:2])
end

function Distributions._rand!(rng::AbstractRNG, s::DynamicalDistribution, x::Vector{<:Array})
    i = rand(rng, 1:length(s.positions))
    x[1] .= select_item(s.velocities, i, s.size)
    x[2] .= select_item(s.positions, i, s.size)
    x
end

select_item(x::Vector, i::Integer, ::NTuple) = x[i]
select_item(x::Sampleable{Univariate}, ::Integer, size::NTuple) = rand(x, size)
select_item(x::Real, ::Integer, ::NTuple) = x
