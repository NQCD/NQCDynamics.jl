using Test
using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.InitialConditions
using Random
using Distributions
using HDF5

a = [rand(Float64, 3, 2) for i=1:10]
d = DynamicalDistribution(a, a, (3,2))
@test eltype(a[1]) == eltype(d)
@test (3,2) == size(d)
@test rand(d) isa Vector{<:Matrix}
@test InitialConditions.maxindex(d) == 10

d = DynamicalDistribution(a, Normal(), (3,2))
@test rand(d) isa Vector{<:Matrix}
@test InitialConditions.maxindex(d) == 10

d = DynamicalDistribution(1, Normal(), (3,2))
@test rand(d) isa Vector{<:Matrix}
@test InitialConditions.maxindex(d) == 1

@testset "BoltzmannVelocityDistribution" begin
    boltz = InitialConditions.BoltzmannVelocityDistribution(1, [1, 2, 3])
    @test length(rand(boltz)) == 3
    @test size(InitialConditions.select_item(boltz, 1, (2, 3))) == (2, 3)
    @test size(InitialConditions.select_item(boltz, 1, (2, 3, 4))) == (2, 3, 4)
end

@testset "HDF5 distribution" begin
    filename = "dat.h5"
    InitialConditions.write_hdf5(filename, a)
    @test InitialConditions.select_item(filename, 1, (3,2)) ∈ a

    d = DynamicalDistribution(filename, 1, (3,2))
    v, r = rand(d)
    @test all(r .== 1)
    @test v ∈ a
    @test InitialConditions.maxindex(d) == 10
    rm(filename)
end
