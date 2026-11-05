using LinearAlgebra
using RingPolymerArrays
using NQCModels
using NQCCalculators

@testset "QuantumModel_Cache" begin
    model = NQCModels.DoubleWell()
    r = rand(1,1)
    cache = NQCCalculators.QuantumModel_Cache(model, 1, r)
    quantities = (:potential, :derivative, :eigen, :adiabatic_derivative, :nonadiabatic_coupling)

    @testset "Potential evaluation" begin
        #check that the position is only changed in cache after get_potential() is called and not when potential() is called
        r = hcat(0.1)
        true_potential = NQCModels.potential(model, r)

        @test cache.potential != true_potential
        #@test getfield(cache, :potential).position != r

        #NQCCalculators.get_potential(cache, r)
        NQCCalculators.update_potential!(cache, r)

        @test cache.potential == true_potential

        # Check that the position for potential is the only one that has been updated.
        @test NQCCalculators.get_derivative(cache, r)[1,1] !== NQCModels.derivative(model, r)
        @test NQCCalculators.get_eigen(cache, r).values !== eigvals(NQCModels.potential(model, r))
        @test abs.(NQCCalculators.get_eigen(cache, r).vectors) !== abs.(eigvecs(NQCModels.potential(model, r)))
        @test NQCCalculators.get_adiabatic_derivative(cache, r)[1] !== cache.eigen.vectors' * NQCModels.derivative(model, r)[1] * cache.eigen.vectors
    end

    @testset "Dependent evaluation" begin
        #Change position value
        r = hcat(0.2)

        #Check that get_:quantities retrieves the wrong values
        @test NQCCalculators.get_potential(cache, r) !== NQCModels.potential(model, r)
        @test NQCCalculators.get_derivative(cache, r)[1] !== NQCModels.derivative(model, r)
        @test NQCCalculators.get_eigen(cache, r).values !== eigvals(NQCModels.potential(model, r))
        @test abs.(NQCCalculators.get_eigen(cache, r).vectors) !== abs.(eigvecs(NQCModels.potential(model, r)))
        @test NQCCalculators.get_adiabatic_derivative(cache, r)[1] !== cache.eigen.vectors' * NQCModels.derivative(model, r)[1] * cache.eigen.vectors

		@debug "Cache info" cache
        #Update entries in the Cache to the new values
        NQCCalculators.update_cache!(cache, r)
		@debug "Cache info" cache

        # Check all the results
        @test cache.potential ≈ NQCModels.potential(model, r)
        @test cache.derivative[1] ≈ NQCModels.derivative(model, r)[1]
        @test cache.eigen.values ≈ eigvals(NQCModels.potential(model, r))
        @test abs.(cache.eigen.vectors) ≈ abs.(eigvecs(NQCModels.potential(model, r)))
        @test isapprox(
            cache.adiabatic_derivative[1],
            cache.eigen.vectors' * NQCModels.derivative(model, r)[1] * cache.eigen.vectors
        )
    end

    allocations = [
        @allocated(NQCCalculators.get_potential(cache, r)),
        @allocated(NQCCalculators.get_derivative(cache, r)),
        @allocated(NQCCalculators.get_eigen(cache, r)),
        @allocated(NQCCalculators.get_adiabatic_derivative(cache, r)),
        @allocated(NQCCalculators.get_nonadiabatic_coupling(cache, r)),
    ]
    if any(allocation > 0 for allocation in allocations)
        @warn "Some get_:quantities functions allocated memory, which should be avoided: $(quantities[findall(allocations .> 0)])"
    end

    @testset "Explict Bath Model" begin
        model = NQCModels.AndersonHolstein(NQCModels.ErpenbeckThoss(;Γ=1.0), TrapezoidalRule(40, -1.0, 1.0))
        r = rand(1,1)
        cache = NQCCalculators.QuantumModel_Cache(model, 1, r)

        NQCCalculators.update_cache!(cache, r)

        @test cache.potential ≈ NQCModels.potential(model, r)
        @test cache.derivative ≈ NQCModels.derivative(model, r)
        @test cache.eigen.values ≈ eigvals(NQCModels.potential(model, r))
        @test abs.(cache.eigen.vectors) ≈ abs.(eigvecs(NQCModels.potential(model, r)))
        @test isapprox(
            cache.adiabatic_derivative[1],
            cache.eigen.vectors' * NQCModels.derivative(model, r)[1] * cache.eigen.vectors
        )
    end
end
