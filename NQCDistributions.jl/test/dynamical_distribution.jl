using NQCDistributions
using Distributions: Normal

@testset "iterate" begin
    d = DynamicalDistribution(Normal(), 2.0, (3,2))
    @test length(collect(d)) == 1

    d = DynamicalDistribution([Normal() for i=1:3, j=1:2], [rand(3,2) for _=1:10], (3,2))
    @test length(collect(d)) == 10
end

@testset "rand" begin
    d = DynamicalDistribution(Normal(), 2.0, (3,2))
    @test rand(d) != rand(d)
    @test rand(d) != rand(d)
end

@testset "getindex" begin
    d = DynamicalDistribution([Normal() for i=1:3, j=1:2], Normal(), (3,2))
    @test d[1] == d[1]
    @test d[10] == d[10]
    @test d[100] == d[100]

    d = DynamicalDistribution([Normal() for i=1:3, j=1:2], [rand(3,2) for _=1:10], (3,2))
    @test d[1] == d[1]
    @test d[10] == d[10]
    @test_throws BoundsError d[100] == d[100]
end
