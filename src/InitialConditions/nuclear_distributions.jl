
export DynamicalDistribution

using Random: AbstractRNG, rand!
using Distributions
using UnitfulAtomic

struct DynamicalVariate <: VariateForm end

struct DynamicalDistribution{V,R,S,N} <: Sampleable{DynamicalVariate,Continuous}
    velocity::V
    position::R
    size::NTuple{S,Int}
    state::N
    type::Symbol
end
DynamicalDistribution(velocity, position, size; state=0, type=:adiabatic) =
    DynamicalDistribution(velocity, position, size, state, type)

Base.eltype(s::DynamicalDistribution{<:Sampleable,R}) where {R} = eltype(s.velocity)
Base.eltype(s::DynamicalDistribution{<:AbstractArray,R} where {R}) = eltype(austrip.(s.velocity[1]))
Base.size(s::DynamicalDistribution) = s.size

function Distributions.rand(rng::AbstractRNG, s::Sampleable{DynamicalVariate})
    Distributions._rand!(rng, s, [Array{eltype(s)}(undef, size(s)) for i=1:2])
end

function Distributions._rand!(rng::AbstractRNG, s::DynamicalDistribution, x::Vector{<:Array})
    i = rand(rng, 1:length(s.position))
    x[1] .= select_item(s.velocity, i, s.size)
    x[2] .= select_item(s.position, i, s.size)
    x
end

pick(s::DynamicalDistribution, i::Integer) = [select_item(s.velocity, i, s.size), select_item(s.position, i, s.size)]

select_item(x::Vector, i::Integer, ::NTuple) = austrip.(x[i])
select_item(x::Sampleable{Univariate}, ::Integer, size::NTuple) = austrip.(rand(x, size))
select_item(x::Real, ::Integer, size::NTuple) = austrip.(fill(x, size))
select_item(x::Matrix, ::Integer, ::NTuple) = austrip.(x)
