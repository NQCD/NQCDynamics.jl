using Test
using NQCDistributions
using NQCDistributions: SampleableComponent
using Distributions: Normal
using Random: Random, SamplerTrivial, Xoshiro
using CUDA

@testset "Extension loads" begin
    @test Base.get_extension(NQCDistributions, :NQCDistributionsCUDAExt) !== nothing
end

if CUDA.functional()
    @testset "rand! into CuArray" begin
        dims = (3, 2)
        rng = Xoshiro(1)

        @testset "UnivariateFill" begin
            a = CUDA.zeros(Float64, dims)
            Random.rand!(rng, a, SamplerTrivial(SampleableComponent(Normal(), dims)))
            @test a isa CuArray
            @test all(!iszero, Array(a))
        end

        @testset "UnivariateArray" begin
            sampleables = [Normal(i, 0.1) for i = 1:dims[1], j = 1:dims[2]]
            a = CUDA.zeros(Float64, dims)
            Random.rand!(rng, a, SamplerTrivial(SampleableComponent(sampleables, dims)))
            @test all(isapprox.(Array(a), [i for i = 1:dims[1], j = 1:dims[2]]; atol = 1))
        end

        @testset "FixedArray" begin
            value = rand(dims...)
            a = CUDA.zeros(Float64, dims)
            Random.rand!(rng, a, SamplerTrivial(SampleableComponent(value, dims)))
            @test Array(a) == value
        end

        @testset "FixedFill" begin
            a = CUDA.zeros(Float64, dims)
            Random.rand!(rng, a, SamplerTrivial(SampleableComponent(2.0, dims)))
            @test all(==(2.0), Array(a))
        end

        @testset "ConfigurationVector" begin
            configurations = [fill(Float64(i), dims) for i = 1:5]
            a = CUDA.zeros(Float64, dims)
            Random.rand!(rng, a, SamplerTrivial(SampleableComponent(configurations, dims)))
            @test Array(a) in configurations
        end
    end
else
    @info "CUDA is not functional, skipping GPU sampling tests."
    @test_skip CUDA.functional()
end
