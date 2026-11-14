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

@testset "frozen_atoms" begin
    # Non-indexable velocity distribution: frozen atoms always have zero velocity
    d = DynamicalDistribution(Normal(), 2.0, (3, 4); frozen_atoms=[1, 3])
    @test d.frozen_atoms == [false, true, false, true]
    s = rand(d)
    @test all(iszero, s.v[:, [1, 3]])
    @test all(!iszero, s.v[:, [2, 4]])
    @test all(iszero, d[5].v[:, [1, 3]])

    # Without frozen atoms the velocities are left untouched
    d = DynamicalDistribution(1.0, 2.0, (3, 4))
    @test all(rand(d).v .== 1.0)

    # Indexable velocity distribution (vector of configurations)
    velocities = [ones(3, 4) for _ in 1:5]
    d = DynamicalDistribution(velocities, 2.0, (3, 4); frozen_atoms=[2])
    @test lastindex(d) == 5
    @test length(d) == 5
    @test length(collect(d)) == 5
    @test all(iszero, d[3].v[:, 2])
    @test all(d[3].v[:, [1, 3, 4]] .== 1.0)
    @test all(iszero, rand(d).v[:, 2])
    @test velocities[3] == ones(3, 4) # Stored configurations are not modified

    # Ring polymer distribution: freezing applies to every bead
    d = DynamicalDistribution(Normal(), 2.0, (3, 4, 5); frozen_atoms=[4])
    s = rand(d)
    @test size(s.v) == (3, 4, 5)
    @test all(iszero, s.v[:, 4, :])
    @test all(!iszero, s.v[:, 1:3, :])
end
