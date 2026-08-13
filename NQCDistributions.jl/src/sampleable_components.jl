
abstract type SampleableComponent end

"""
    UnivariateFill{S<:Sampleable{Univariate}}
  
Fill all degrees of freedom from single Univariate distribution.
"""
struct UnivariateFill{S<:Sampleable{Univariate}} <: SampleableComponent
    sampleable::S
    dims::Dims{2}
end
function Random.rand(rng::AbstractRNG, d::SamplerTrivial{<:UnivariateFill})
    return rand(rng, d[].sampleable, d[].dims)
end
function Random.rand!(rng::AbstractRNG, a::AbstractMatrix, d::SamplerTrivial{<:UnivariateFill})
    return Random.rand!(rng, d[].sampleable, a)
end
Base.eltype(d::UnivariateFill) = Matrix{eltype(d.sampleable)}
Base.size(d::UnivariateFill) = d.dims
isindexable(::UnivariateFill) = false

"""
    UnivariateArray{N,S<:Sampleable{Univariate}}

Fill each degree of freedom from a different Univariate distribution.

The size of the matrix of sampleables should match the system size.
"""
struct UnivariateArray{N,S<:Sampleable{Univariate}} <: SampleableComponent
    sampleable::Array{S,N}
end
function Random.rand(rng::AbstractRNG, d::SamplerTrivial{<:UnivariateArray})
    out = similar(d[].sampleable, Union{unique(eltype.(d[].sampleable))...})
    return Random.rand!(rng, out, d)
end
function Random.rand!(rng::AbstractRNG, a::AbstractArray, d::SamplerTrivial{<:UnivariateArray})
    for I in eachindex(a, d[].sampleable)
        a[I] = rand(rng, d[].sampleable[I])
    end
    return a
end
Base.eltype(d::UnivariateArray{N}) where {N} = Array{eltype(eltype(d.sampleable)),N}
Base.size(d::UnivariateArray) = size(d.sampleable)
isindexable(::UnivariateArray) = false

"""
    FixedArray{S<:AbstractArray}

Return the same configuration every time.
"""
struct FixedArray{S<:AbstractArray} <: SampleableComponent
    value::S
end
function Random.rand(::AbstractRNG, d::SamplerTrivial{<:FixedArray})
    return copy(d[].value)
end
function Random.rand!(::AbstractRNG, a::AbstractArray, d::SamplerTrivial{<:FixedArray})
    copy!(a, d[].value)
    return a
end
Base.eltype(d::FixedArray) = typeof(d.value)
Base.size(d::FixedArray) = size(d.value)
isindexable(::FixedArray) = false

"""
    FixedFill{S<:Real}

Fill all degrees of freedom with the same value every time.
"""
struct FixedFill{S<:Real} <: SampleableComponent
    value::S
    dims::Dims{2}
end
function Random.rand(::AbstractRNG, d::SamplerTrivial{<:FixedFill})
    return fill(d[].value, d[].dims)
end
function Random.rand!(::AbstractRNG, a::AbstractMatrix, d::SamplerTrivial{<:FixedFill})
    fill!(a, d[].value)
    return a
end
Base.eltype(d::FixedFill) = Matrix{typeof(d.value)}
Base.size(d::FixedFill) = d.dims
isindexable(::FixedFill) = false

"""
    ConfigurationVector{S<:AbstractVector}

Sample from a provided vector of configurations.
"""
struct ConfigurationVector{S<:AbstractVector} <: SampleableComponent
    configurations::S
end
function Random.rand(rng::AbstractRNG, d::SamplerTrivial{<:ConfigurationVector})
    return copy(rand(rng, d[].configurations))
end
function Random.rand!(rng::AbstractRNG, a::AbstractArray, d::SamplerTrivial{<:ConfigurationVector})
    copy!(a, rand(rng, d[].configurations))
    return a
end
Base.eltype(d::ConfigurationVector) = eltype(d.configurations)
Base.size(d::ConfigurationVector) = size(first(d.configurations))
Base.getindex(d::ConfigurationVector, i) = copy(d.configurations[i])
Base.firstindex(d::ConfigurationVector) = firstindex(d.configurations)
Base.lastindex(d::ConfigurationVector) = lastindex(d.configurations)
isindexable(::ConfigurationVector) = true

"""
    RingPolymerWrapper{S}

Wrap other distributions to convert them to ring polymer distributions.
"""
struct RingPolymerWrapper{S} <: SampleableComponent
    sampleable::S
    dims::Dims{3}
    classical::Vector{Int}
end
function RingPolymerWrapper(sampleable, nbeads, classical)
    return RingPolymerWrapper(sampleable, (size(sampleable)..., nbeads), classical)
end
function Random.rand(rng::AbstractRNG, d::SamplerTrivial{<:RingPolymerWrapper})
    out = RingPolymerArray{eltype(eltype(d[].sampleable))}(undef, d[].dims; d[].classical)
    Random.rand!(rng, out, d)
    return out
end
function Random.rand!(rng::AbstractRNG, a::RingPolymerArray, d::SamplerTrivial{<:RingPolymerWrapper})
    for bead in eachbead(a)
        Random.rand!(rng, bead, d[].sampleable)
    end
    return a
end
Base.eltype(d::RingPolymerWrapper) = RingPolymerArray{eltype(eltype(d.sampleable))}
Base.size(d::RingPolymerWrapper) = d.dims
function Base.getindex(d::RingPolymerWrapper, i)
    out = RingPolymerArray{eltype(eltype(d.sampleable))}(undef, d.dims; d.classical)
    for bead in eachbead(out)
        copy!(bead, d.sampleable[i])
    end
    return out
end
Base.firstindex(d::RingPolymerWrapper) = firstindex(d.sampleable)
Base.lastindex(d::RingPolymerWrapper) = lastindex(d.sampleable)
isindexable(d::RingPolymerWrapper) = isindexable(d.sampleable)

"""
    SampleableComponent(sampleable, dims)

Converts a general `sampleable` that provides configurations into one of the component types defined below. 
`dims` should be the size of the desired samples and must be consistent with the provided sampleable.
"""
function SampleableComponent end

function SampleableComponent(sampleable::SampleableComponent, dims::Dims{3}, classical)
    checkdims(size(sampleable), (dims[1], dims[2]))
    return RingPolymerWrapper(sampleable, dims[3], classical)
end

function SampleableComponent(sampleable::SampleableComponent, dims::Dims{2})
    checkdims(size(sampleable), dims)
    return sampleable
end

function SampleableComponent(sampleable::Sampleable{Univariate}, dims::Dims{2})
    return UnivariateFill(sampleable, dims)
end
function SampleableComponent(sampleable::Sampleable{Univariate}, dims::Dims{3}, classical)
    return RingPolymerWrapper(UnivariateFill(sampleable, (dims[1], dims[2])), dims[3], classical)
end

function SampleableComponent(sampleable::AbstractArray{<:Sampleable{Univariate}}, dims::Dims{2})
    checkdims(size(sampleable), dims)
    return UnivariateArray(sampleable)
end
function SampleableComponent(sampleable::AbstractArray{<:Sampleable{Univariate}}, dims::Dims{3}, classical)
    checkdims(size(sampleable), (dims[1], dims[2]))
    return RingPolymerWrapper(UnivariateArray(sampleable), dims[3], classical)
end

function SampleableComponent(sampleable::AbstractArray{<:Real}, dims::Dims)
    checkdims(size(sampleable), dims)
    return FixedArray(sampleable)
end

function SampleableComponent(sampleable::Real, dims::Dims{2})
    return FixedFill(sampleable, dims)
end
function SampleableComponent(sampleable::Real, dims::Dims{3}, classical)
    return RingPolymerWrapper(FixedFill(sampleable, (dims[1], dims[2])), dims[3], classical)
end

function SampleableComponent(sampleable::AbstractVector{<:AbstractArray}, dims::Dims)
    checkdims(size(sampleable[1]), dims)
    return ConfigurationVector(sampleable)
end
function SampleableComponent(sampleable::AbstractVector{<:AbstractArray}, dims::Dims{3}, classical)
    checkdims(size(sampleable[1]), dims)
    return ConfigurationVector(sampleable)
end

function SampleableComponent(sampleable::RingPolymerWrapper, dims::Dims{3}, classical)
    checkdims(size(sampleable), dims)
    return sampleable
end

checkdims(sz, dims) = sz == dims || throw(DimensionMismatch("Size of sampleable does not match provided dims: $sz != $dims."))
