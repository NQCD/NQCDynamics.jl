using Test
using LinearAlgebra
using RingPolymerArrays
using NQCModels
using NQCCalculators

# Allocation-based tests aren't great - change these values if there's a good reason to.
# Eigenvalue allocations limit: Regardless of maths backend and Julia version, eigenvalues should take fewer allocations than this or something might be wrong
const MAX_EIGEN_ALLOCATIONS = 60000

@testset "General constructors" begin #reasonable basic tests to ensure that multiple dispatch is correctly working to construct the correct types of Calculators
    model = NQCModels.DoubleWell()
    @test NQCCalculators.Create_Cache(model, 1, hcat(1.0)) isa NQCCalculators.QuantumModel_Cache
    @test NQCCalculators.Create_Cache(model, 1, 2, [1.0;;;]) isa NQCCalculators.RingPolymer_QuantumModel_Cache

    model = NQCModels.Free()
    @test NQCCalculators.Create_Cache(model, 1, hcat(1.0)) isa NQCCalculators.ClassicalModel_Cache
    @test NQCCalculators.Create_Cache(model, 1, 2, [1.0;;;]) isa NQCCalculators.RingPolymer_ClassicalModel_Cache

    model = NQCModels.AndersonHolstein(MiaoSubotnik(;Γ=1.0), TrapezoidalRule(10, -1.0, 1.0))
    @test NQCCalculators.Create_Cache(model, 1, hcat(1.0)) isa NQCCalculators.QuantumModel_Cache

    model = NQCModels.CompositeFrictionModel(Free(), RandomFriction(1))
    @test NQCCalculators.Create_Cache(model, 1, hcat(1.0)) isa NQCCalculators.ClassicalFrictionModel_Cache
    @test NQCCalculators.Create_Cache(model, 1, 2, [1.0;;;]) isa NQCCalculators.RingPolymer_ClassicalFrictionModel_Cache
end
