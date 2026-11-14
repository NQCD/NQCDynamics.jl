using NQCDistributions
using Distributions: Normal
using Random: rand!
using RingPolymerArrays: RingPolymerArray

using NQCDistributions:
    SampleableComponent,
    UnivariateArray,
    UnivariateFill,
    FixedArray,
    FixedFill,
    ConfigurationVector,
    RingPolymerWrapper

@testset "UnivariateFill" begin
    d = UnivariateFill(Normal(), (3,2))
    @test eltype(d) == Matrix{Float64}
    @test rand(d) isa AbstractMatrix
    out = zeros(3,2)
    rand!(out, d)
    @test all(out .!== 0.0)
    @test d == SampleableComponent(Normal(), (3,2))
end

@testset "UnivariateArray" begin
    d = UnivariateArray([Normal() for i=1:3, j=1:2])
    @test eltype(d) == Matrix{Float64}
    @test rand(d) isa AbstractMatrix
    out = zeros(3,2)
    rand!(out, d)
    @test all(out .!== 0.0)
    @test SampleableComponent([Normal() for i=1:3, j=1:2], (3,2)) isa UnivariateArray
    @test_throws DimensionMismatch SampleableComponent([Normal() for i=1:3, j=1:2], (2,2))

    d = UnivariateArray([Normal() for i=1:3, j=1:2, k=1:4])
    @test eltype(d) == Array{Float64,3}
    @test rand(d) isa AbstractArray
    out = zeros(3,2,4)
    rand!(out, d)
    @test all(out .!== 0.0)
    @test SampleableComponent([Normal() for i=1:3, j=1:2], (3,2,4), Int[]) isa RingPolymerWrapper
end

@testset "FixedArray" begin
    d = FixedArray([1 2; 3 4; 5 6])
    @test eltype(d) == Matrix{Int}
    @test rand(d) isa AbstractMatrix
    out = zeros(3,2)
    rand!(out, d)
    @test out == [1.0 2.0; 3.0 4.0; 5.0 6.0]
    @test SampleableComponent([1 2; 3 4; 5 6], (3,2)) isa FixedArray
    @test_throws DimensionMismatch SampleableComponent([1 2; 3 4; 5 6], (3,3))

    a = rand(3,2,3)
    d = FixedArray(a)
    @test eltype(d) == Array{Float64,3}
    @test rand(d) isa AbstractArray{Float64,3}
    out = zeros(3,2,3)
    rand!(out, d)
    @test out == a
    @test SampleableComponent(rand(3,2,3), (3,2,3)) isa FixedArray
    @test_throws DimensionMismatch SampleableComponent(rand(3,2,3), (3,3,4))
end

@testset "FixedFill" begin
    d = FixedFill(1.0, (3,2))
    @test eltype(d) == Matrix{Float64}
    @test rand(d) isa AbstractMatrix
    out = zeros(3,2)
    rand!(out, d)
    @test out == ones(3,2)
    @test SampleableComponent(1.0, (3,2)) isa FixedFill
    @test SampleableComponent(1.0, (3,2,3), Int[]) isa RingPolymerWrapper
end

@testset "ConfigurationVector" begin
    d = ConfigurationVector([ones(3,2) for _=1:3])
    @test eltype(d) == Matrix{Float64}
    @test rand(d) isa AbstractMatrix
    @test d[1] isa AbstractMatrix
    out = zeros(3,2)
    rand!(out, d)
    @test out == ones(3,2)
    @test SampleableComponent([ones(3,2) for _=1:3], (3,2)) isa ConfigurationVector
    @test_throws DimensionMismatch SampleableComponent([ones(3,2) for _=1:3], (3,2,3))

    d = ConfigurationVector([ones(3,2,3) for _=1:3])
    @test eltype(d) == Array{Float64,3}
    @test rand(d) isa AbstractArray{Float64,3}
    out = zeros(3,2,3)
    rand!(out, d)
    @test out == ones(3,2,3)
    @test SampleableComponent([ones(3,2,3) for _=1:3], (3,2,3)) isa ConfigurationVector
    @test_throws DimensionMismatch SampleableComponent([ones(3,2,3) for _=1:3], (3,2,2))
end

@testset "RingPolymerWrapper" begin
    nbeads = 4
    classical = [1]

    @testset "UnivariateFill" begin
        dsingle = UnivariateFill(Normal(), (3,2))
        d = RingPolymerWrapper(dsingle, nbeads, classical)
        @test eltype(d) <: RingPolymerArray
        @test rand(d) isa RingPolymerArray
    end

    @testset "UnivariateArray" begin
        dsingle = UnivariateArray([Normal() for i=1:3, j=1:2])
        d = RingPolymerWrapper(dsingle, nbeads, classical)
        @test eltype(d) <: RingPolymerArray
        @test rand(d) isa RingPolymerArray
    end

    @testset "FixedArray" begin
        dsingle = UnivariateArray([Normal() for i=1:3, j=1:2])
        d = RingPolymerWrapper(dsingle, nbeads, classical)
        @test eltype(d) <: RingPolymerArray
        @test rand(d) isa RingPolymerArray
    end

    @testset "FixedFill" begin
        dsingle = FixedFill(1.0, (3,2))
        d = RingPolymerWrapper(dsingle, nbeads, classical)
        @test eltype(d) <: RingPolymerArray
        @test rand(d) isa RingPolymerArray
        @test rand(d) == ones(3,2,4)
    end

    @testset "ConfigurationVector" begin
        dsingle = ConfigurationVector([ones(3,2) for _=1:3])
        d = RingPolymerWrapper(dsingle, nbeads, classical)
        @test eltype(d) <: RingPolymerArray
        @test d[1] isa RingPolymerArray
        @test rand(d) isa RingPolymerArray
        @test rand(d) == ones(3,2,4)
    end
end

