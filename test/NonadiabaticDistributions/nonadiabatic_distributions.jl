
using Test
using NonadiabaticMolecularDynamics
using Distributions: Normal

@testset "ElectronicDistribution" begin
    diabatic = NonadiabaticDistributions.Diabatic()

    @testset "SingleState" begin
        @test NonadiabaticDistributions.SingleState(1) === NonadiabaticDistributions.SingleState(1, diabatic)
    end

    @testset "ElectronicPopulation" begin
        dist = NonadiabaticDistributions.ElectronicPopulation([1, 0])
        dist2 = NonadiabaticDistributions.ElectronicPopulation([1, 0], diabatic)
        @test dist.statetype === dist2.statetype
        @test dist.populations == dist2.populations
    end
end

@testset "CombinedDistribution" begin
    electronic = NonadiabaticDistributions.SingleState(1)
    nuclear = DynamicalDistribution(Normal(), Normal(), (3, 2))
    combined = electronic * nuclear
    @test combined isa NonadiabaticDistributions.NonadiabaticDistribution
end
