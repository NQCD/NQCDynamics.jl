using LinearAlgebra
using RingPolymerArrays
using NQCModels
using NQCCalculators

@testset "ClassicalModel_Cache" begin
    model = NQCModels.Harmonic()
    r = rand(1, 5)
    cache = NQCCalculators.ClassicalModel_Cache(model, 5, r)

    #Check that get_:quantities retrieves the wrong values
    @test NQCCalculators.get_potential(cache, r)[1] !== NQCModels.potential(model, r)
    @test NQCCalculators.get_derivative(cache, r) !== NQCModels.derivative(model, r)
	@debug "Cache info" cache
    NQCCalculators.update_cache!(cache, r)
	@debug "Cache info" cache
    # Check all the results

    @test cache.potential[1] ≈ NQCModels.potential(model, r)
    @test cache.derivative ≈ NQCModels.derivative(model, r)
end
