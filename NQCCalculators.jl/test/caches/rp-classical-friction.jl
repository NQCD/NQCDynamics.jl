using Test
using LinearAlgebra
using RingPolymerArrays
using NQCModels
using NQCCalculators

# Allocation-based tests aren't great - change these values if there's a good reason to.
# Eigenvalue allocations limit: Regardless of maths backend and Julia version, eigenvalues should take fewer allocations than this or something might be wrong
const MAX_EIGEN_ALLOCATIONS = 60000

@testset "RingPolymer_ClassicalFrictionModel_Cache" begin
    model = CompositeFrictionModel(Free(), RandomFriction(1))
    r = rand(1,1,10)
    cache = NQCCalculators.RingPolymer_ClassicalFrictionModel_Cache(model, 1, 10, r)
	@debug "Cache info" cache
    NQCCalculators.get_potential(cache, r)
    NQCCalculators.get_derivative(cache, r)
    NQCCalculators.get_friction(cache, r)

    NQCCalculators.evaluate_potential(cache, r)
    NQCCalculators.evaluate_derivative(cache, r)
    NQCCalculators.evaluate_friction(cache, r)
	@debug "Cache info" cache
end
