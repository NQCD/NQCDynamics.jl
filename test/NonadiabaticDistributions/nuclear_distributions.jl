using Test
using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.InitialConditions
using Random
using Distributions
using HDF5
using ComponentArrays

a = [rand(Float64, 3, 2) for i=1:100]
d = DynamicalDistribution(a, a, (3,2))
@test (3,2) == size(d)
@test rand(d) isa ComponentVector
@test NonadiabaticDistributions.maxindex(d) == 100

d = DynamicalDistribution(a, Normal(), (3,2))
@test rand(d) isa ComponentVector
@test NonadiabaticDistributions.maxindex(d) == 100

d = DynamicalDistribution(1, Normal(), (3,2))
@test rand(d) isa ComponentVector
@test NonadiabaticDistributions.maxindex(d) == 1

@testset "BoltzmannVelocityDistribution" begin
    boltz = NonadiabaticDistributions.BoltzmannVelocityDistribution(1, [1, 2, 3], (2,3))
    @test size(rand(boltz)) === boltz.size
    boltz = NonadiabaticDistributions.BoltzmannVelocityDistribution(1, [1, 2, 3], (2,3,3))
    @test size(rand(boltz)) === boltz.size
end

@testset "HDF5 distribution" begin
    filename = "dat.h5"
    NonadiabaticDistributions.write_hdf5(filename, a)
    @test NonadiabaticDistributions.select_item(filename, 1, (3,2)) ∈ a

    d = DynamicalDistribution(filename, 1, (3,2))
    u = rand(d)
    @test all(u.r .== 1)
    @test u.v ∈ a
    @test NonadiabaticDistributions.maxindex(d) == 100
    rm(filename)
end
