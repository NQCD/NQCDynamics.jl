
export PositionDistribution
export PhasespaceDistribution

using Random: AbstractRNG
using Distributions

struct PositionDistribution{T} <: Sampleable{Matrixvariate,Continuous}
    positions::Vector{Matrix{T}}
end

struct PhasespaceVariate <: VariateForm end

struct PhasespaceDistribution{T} <: Sampleable{PhasespaceVariate,Continuous}
    positions::Vector{Matrix{T}}
    momenta::Vector{Matrix{T}}
end

NuclearDistributions{T} = Union{PositionDistribution{T}, PhasespaceDistribution{T}}

Base.eltype(::NuclearDistributions{T}) where {T} = T
Base.size(s::NuclearDistributions{T}) where {T} = size(s.positions[1])

function Distributions._rand!(rng::AbstractRNG, s::PositionDistribution{T}, x::DenseMatrix{T}) where T<:AbstractFloat
    x .= rand(rng, s.positions)
end

function Distributions.rand(rng::AbstractRNG, s::Sampleable{PhasespaceVariate})
    Distributions._rand!(rng, s, Phasespace([Matrix{eltype(s)}(undef, size(s)) for i=1:2]...))
end

function Distributions._rand!(rng::AbstractRNG, s::PhasespaceDistribution{T}, x::Phasespace{T}) where T<:AbstractFloat
    i = rand(rng, 1:length(s.positions))
    get_positions(x) .= s.positions[i]
    get_momenta(x) .= s.momenta[i]
    x
end
