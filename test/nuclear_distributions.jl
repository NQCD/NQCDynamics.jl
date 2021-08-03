using Test
using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.InitialConditions
using Random
using Distributions

a = [rand(Float64, 3, 2) for i=1:10]
d = DynamicalDistribution(a, a, (3,2))
@test eltype(a[1]) == eltype(d)
@test (3,2) == size(d)
@test rand(d) isa Vector{<:Matrix}

d = DynamicalDistribution(a, Normal(), (3,2))
@test rand(d) isa Vector{<:Matrix}

d = DynamicalDistribution(1, Normal(), (3,2))
@test rand(d) isa Vector{<:Matrix}

@testset "BoltzmannVelocityDistribution" begin
    boltz = InitialConditions.BoltzmannVelocityDistribution(1, [1, 2, 3])
    @test length(rand(boltz)) == 3
    @test size(InitialConditions.select_item(boltz, 1, (2, 3))) == (2, 3)
    @test size(InitialConditions.select_item(boltz, 1, (2, 3, 4))) == (2, 3, 4)
end
