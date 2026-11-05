using LinearAlgebra
using RingPolymerArrays
using NQCModels
using NQCCalculators

@testset "RingPolymer_QuantumModel_Cache" begin
    model = NQCModels.DoubleWell()
    r = rand(1,1,10)
    reset_cache() = NQCCalculators.RingPolymer_QuantumModel_Cache(model, 1, 10, r)

    r_centroid = rand(1,1)
    standard_quantities = (:potential, :derivative, :eigen, :adiabatic_derivative, :nonadiabatic_coupling)
    centroid_quantities = (:centroid, :centroid_potential, :centroid_derivative,
        :centroid_eigen, :centroid_adiabatic_derivative, :centroid_nonadiabatic_coupling)

    @testset "Potential evaluation" begin
        cache = reset_cache()
        true_potential = [NQCModels.potential(model, r[:,:,i]) for i=1:10]

        @test cache.potential != true_potential
        #@test getfield(cache, :potential).position != r

        NQCCalculators.evaluate_potential(cache, r)

        @test cache.potential != true_potential

        NQCCalculators.update_potential!(cache, r)

        @test cache.potential == true_potential

        # Check that the potential is the only one that has been updated.
        for i in NQCCalculators.beads(cache)
            @test NQCCalculators.get_derivative(cache, r)[i] !== NQCModels.derivative(model, r[:,:,i])[1]
            @test NQCCalculators.get_eigen(cache, r)[i].values !== eigvals(NQCModels.potential(model, r[:,:,i]))
            @test abs.(NQCCalculators.get_eigen(cache, r)[i].vectors) !== abs.(eigvecs(NQCModels.potential(model, r[:,:,i])))
            @test NQCCalculators.get_adiabatic_derivative(cache, r)[i] !== cache.eigen[i].vectors' * NQCModels.derivative(model, r[:,:,i])[1] * cache.eigen[i].vectors
        end
		@debug "Cache info" cache
        NQCCalculators.update_cache!(cache, r)
		@debug "Cache info" cache
        # Check that everything has been updated
        for i in NQCCalculators.beads(cache)
            @test NQCCalculators.get_derivative(cache, r)[i] == NQCModels.derivative(model, r[:,:,i])[1]
            @test NQCCalculators.get_eigen(cache, r)[i].values == eigvals(NQCModels.potential(model, r[:,:,i]))
            @test abs.(NQCCalculators.get_eigen(cache, r)[i].vectors) == abs.(eigvecs(NQCModels.potential(model, r[:,:,i])))
            @test NQCCalculators.get_adiabatic_derivative(cache, r)[i] == cache.eigen[i].vectors' * NQCModels.derivative(model, r[:,:,i])[1] * cache.eigen[i].vectors
        end
    end

    @testset "Dependent evaluation" begin
        cache = reset_cache()

        NQCCalculators.update_cache!(cache, r)

        # Check all the results
        @test cache.potential ≈ [NQCModels.potential(model, r[:,:,i]) for i=1:10]
        @test cache.derivative ≈ [NQCModels.derivative(model, hcat(r[k,j,i]))[1] for k=1:1, j=1:1, i=1:10]
        @test [cache.eigen[i].values for i=1:10] ≈ eigvals.(cache.potential)
        @test [abs.(cache.eigen[i].vectors) for i=1:10] ≈ [abs.(eigvecs(cache.potential[i])) for i=1:10]
        @test isapprox(
            cache.adiabatic_derivative[1],
            cache.eigen[1].vectors' * NQCModels.derivative(model, r[:,:,1])[1] * cache.eigen[1].vectors
        )
    end

    @testset "Centroid dependent evaluation" begin
        cache = reset_cache()
        r_centroid = RingPolymerArrays.get_centroid(r)


        NQCCalculators.update_cache!(cache, r)

        # Check all the results
        @test cache.centroid_potential ≈ NQCModels.potential(model, r_centroid)
        @test cache.centroid_derivative[1] ≈ NQCModels.derivative(model, r_centroid)[1]
        @test cache.centroid_eigen.values ≈ eigvals(NQCModels.potential(model, r_centroid))
        @test abs.(cache.centroid_eigen.vectors) ≈ abs.(eigvecs(NQCModels.potential(model, r_centroid)))
        @test isapprox(
            cache.centroid_adiabatic_derivative[1],
            cache.centroid_eigen.vectors' * NQCModels.derivative(model, r_centroid)[1] * cache.centroid_eigen.vectors
        )
    end

    @testset "Extra quantities" begin
        cache = reset_cache()

        NQCCalculators.update_cache!(cache, r)
        NQCCalculators.get_traceless_derivative(cache, r)
        NQCCalculators.get_traceless_adiabatic_derivative(cache, r)

        @test cache.V̄[1] ≈ tr(cache.potential[1]) / nstates(model)
        @test tr(cache.traceless_potential[1]) ≈ 0 atol=eps(Float64)
        @test cache.D̄[1] ≈ tr(cache.derivative[1]) / nstates(model)
        @test tr(cache.traceless_derivative[1]) ≈ 0 atol=eps(Float64)
        @test cache.traceless_adiabatic_derivative[1] ≈ cache.adiabatic_derivative[1] .- Diagonal(fill(cache.D̄[1], nstates(model)))
    end

    @testset "Zero allocations" begin

        cache = reset_cache()

        # Call all functions to ensure they have been compiled
        NQCCalculators.evaluate_potential(cache, r)
        NQCCalculators.evaluate_derivative(cache, r)
        NQCCalculators.evaluate_eigen(cache, r)
        NQCCalculators.evaluate_adiabatic_derivative(cache, r)
        NQCCalculators.evaluate_nonadiabatic_coupling(cache, r)

        NQCCalculators.evaluate_traceless_potential(cache, r)
        NQCCalculators.evaluate_V̄(cache, r)
        NQCCalculators.evaluate_traceless_derivative(cache, r)
        NQCCalculators.evaluate_D̄(cache, r)
        NQCCalculators.evaluate_traceless_adiabatic_derivative(cache, r)

        NQCCalculators.evaluate_centroid(cache, r)
        NQCCalculators.evaluate_centroid_potential(cache, r)
        NQCCalculators.evaluate_centroid_derivative(cache, r)
        NQCCalculators.evaluate_centroid_eigen(cache, r)
        NQCCalculators.evaluate_centroid_adiabatic_derivative(cache, r)
        NQCCalculators.evaluate_centroid_nonadiabatic_coupling(cache, r)

        cache = reset_cache()
    end


    @testset "RingPolymer explict bath model" begin
        model = NQCModels.AndersonHolstein(MiaoSubotnik(;Γ=1.0), TrapezoidalRule(40, -1.0, 1.0))
        cache = NQCCalculators.RingPolymer_QuantumModel_Cache(model, 1, 10, hcat(1.0))
        r = rand(1, 1, 10)
		@debug "Cache info" cache
        NQCCalculators.update_cache!(cache, r)
		@debug "Cache info" cache
        NQCCalculators.evaluate_potential(cache, r)
        NQCCalculators.evaluate_derivative(cache, r)
        NQCCalculators.evaluate_eigen(cache, r)
        NQCCalculators.evaluate_adiabatic_derivative(cache, r)
        NQCCalculators.evaluate_nonadiabatic_coupling(cache, r)
    end
end
