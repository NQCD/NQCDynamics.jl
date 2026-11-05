using LinearAlgebra
using RingPolymerArrays
using NQCModels
using NQCCalculators

@testset "RingPolymer_ClassicalModel_Cache" begin
    model = NQCModels.Harmonic()
    r = rand(1, 5, 10)
    cache = NQCCalculators.RingPolymer_ClassicalModel_Cache(model, 5, 10, r)
    centroid = RingPolymerArrays.get_centroid(r)

    #Check that get_:quantities retrieves the wrong values
    for I in NQCCalculators.beads(cache)
        @test NQCCalculators.get_potential(cache, r)[I] !== hcat(NQCModels.potential(model, r[:,:,I]))
        @test NQCCalculators.get_derivative(cache, r)[:,:,I] !== NQCModels.derivative(model, r[:,:,I])
    end
    @test NQCCalculators.get_centroid(cache, r) !== centroid
    @test NQCCalculators.get_centroid_potential(cache, r) !== hcat(NQCModels.potential(model, centroid))
    @test NQCCalculators.get_centroid_derivative(cache, r) !== NQCModels.derivative(model, centroid)

	@debug "Cache info" cache
    #Update entries in the Cache to the values pertaining to current position
    NQCCalculators.update_cache!(cache, r)
	@debug "Cache info" cache

    # Check all the results
    for I in NQCCalculators.beads(cache)
        @test cache.potential[I] ≈ hcat(NQCModels.potential(model, r[:,:,I]))
        @test cache.derivative[:,:,I] ≈ NQCModels.derivative(model, r[:,:,I])
    end
    @test cache.centroid ≈ centroid
    @test cache.centroid_potential ≈ hcat(NQCModels.potential(model, centroid))
    @test cache.centroid_derivative ≈ NQCModels.derivative(model, centroid)
end
