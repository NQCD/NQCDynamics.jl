module NQCDistributionsCUDAExt

using NQCDistributions
using NQCDistributions: SampleableComponent, UnivariateFill, UnivariateArray, FixedArray, FixedFill, ConfigurationVector
using CUDA
import Random
using Random: AbstractRNG, SamplerTrivial

# When sampling into a CUDA array, perform the sampling on the CPU first,
# then transfer the result to the GPU. This avoids scalar indexing on CUDA arrays.

function Random.rand!(rng::AbstractRNG, a::AnyCuArray, d::SamplerTrivial{<:UnivariateFill})
    cpu_array = rand(rng, d[].sampleable, d[].dims)
    copyto!(a, cpu_array)
    return a
end

function Random.rand!(rng::AbstractRNG, a::AnyCuArray, d::SamplerTrivial{<:UnivariateArray})
    cpu_array = similar(Array(a))
    for I in eachindex(cpu_array, d[].sampleable)
        cpu_array[I] = rand(rng, d[].sampleable[I])
    end
    copyto!(a, cpu_array)
    return a
end

function Random.rand!(rng::AbstractRNG, a::AnyCuArray, d::SamplerTrivial{<:FixedArray})
    copyto!(a, d[].value)
    return a
end

function Random.rand!(rng::AbstractRNG, a::AnyCuArray, d::SamplerTrivial{<:FixedFill})
    fill!(a, d[].value)
    return a
end

function Random.rand!(rng::AbstractRNG, a::AnyCuArray, d::SamplerTrivial{<:ConfigurationVector})
    copyto!(a, rand(rng, d[].configurations))
    return a
end

end
