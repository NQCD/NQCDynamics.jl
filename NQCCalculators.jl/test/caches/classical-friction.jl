using LinearAlgebra
using RingPolymerArrays
using NQCModels
using NQCCalculators

@testset "ClassicalFrictionModel_Cache" begin
    model = CompositeFrictionModel(Free(), RandomFriction(1))
    r = rand(1,1)
    cache = NQCCalculators.ClassicalFrictionModel_Cache(model, 1, r)
	@debug "Cache info" cache
    NQCCalculators.get_potential(cache, r)
    NQCCalculators.get_derivative(cache, r)
    NQCCalculators.get_friction(cache, r)

    NQCCalculators.evaluate_potential(cache, r)
    NQCCalculators.evaluate_derivative(cache, r)
    NQCCalculators.evaluate_friction(cache, r)
	@debug "Cache info" cache
end
