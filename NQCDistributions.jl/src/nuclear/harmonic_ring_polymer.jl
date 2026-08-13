
struct PositionHarmonicRingPolymer{T} <: SampleableComponent
    normal_mode_distribution::UnivariateArray{3,Normal{T}}
    transformation::NormalModeTransformation{T}
    classical::Vector{Int}
end

"""
    PositionHarmonicRingPolymer{T}(ω, β, m, dims::Dims{3}; centre=0, classical=Int[])

Ring polymer position distribution in a 1D harmonic potential
"""
function PositionHarmonicRingPolymer{T}(ω, β, m, dims::Dims{3}; centre=0, classical=Int[]) where {T}

    nbeads = dims[3]

    βₙ = β / nbeads
    ωₙ = 1 / βₙ
    ωₖ = 2ωₙ*sin.((0:nbeads-1)*π/nbeads)
    σ⁻¹ = sqrt.(βₙ*m .* (ω^2 .+ ωₖ.^2))
    σ = T.(1 ./ σ⁻¹)

    μ = T[centre * sqrt(nbeads); zeros(nbeads-1)]

    normal_mode_distribution = UnivariateArray([Normal{T}(μ[I[3]], σ[I[3]]) for I in CartesianIndices(dims)])
    normal_mode_transformation = NormalModeTransformation{T}(nbeads)
    
    PositionHarmonicRingPolymer{T}(normal_mode_distribution, normal_mode_transformation, classical)
end

function Random.rand(rng::AbstractRNG, d::SamplerTrivial{PositionHarmonicRingPolymer{T}}) where {T}
    a = RingPolymerArray{T}(undef, size(d[].normal_mode_distribution); classical=d[].classical)
    Random.rand!(rng, a, d)
    return a
end
function Random.rand!(rng::AbstractRNG, a::AbstractArray, d::SamplerTrivial{<:PositionHarmonicRingPolymer})
    Random.rand!(rng, a, d[].normal_mode_distribution)
    transform_from_normal_modes!(a, d[].transformation)
    return a
end
Base.eltype(::PositionHarmonicRingPolymer{T}) where {T} = RingPolymerArray{T}
Base.size(d::PositionHarmonicRingPolymer) = size(d.normal_mode_distribution)
isindexable(::PositionHarmonicRingPolymer) = false

function SampleableComponent(sampleable::PositionHarmonicRingPolymer, dims::Dims{3}, classical)
    checkdims(size(sampleable), dims)
    sampleable.classical == classical || throw(error("`classical` vectors do not match."))
    return sampleable
end
