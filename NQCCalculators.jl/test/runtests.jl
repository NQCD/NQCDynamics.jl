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
    @test NQCCalculators.Create_Cache(model, 1, 2, hcat(1.0)) isa NQCCalculators.RingPolymer_QuantumModel_Cache

    model = NQCModels.Free()
    @test NQCCalculators.Create_Cache(model, 1, hcat(1.0)) isa NQCCalculators.ClassicalModel_Cache
    @test NQCCalculators.Create_Cache(model, 1, 2, hcat(1.0)) isa NQCCalculators.RingPolymer_ClassicalModel_Cache

    model = NQCModels.AndersonHolstein(MiaoSubotnik(;Γ=1.0), TrapezoidalRule(10, -1.0, 1.0))
    @test NQCCalculators.Create_Cache(model, 1, hcat(1.0)) isa NQCCalculators.QuantumModel_Cache

    model = NQCModels.CompositeFrictionModel(Free(), RandomFriction(1))
    @test NQCCalculators.Create_Cache(model, 1, hcat(1.0)) isa NQCCalculators.ClassicalFrictionModel_Cache
    @test NQCCalculators.Create_Cache(model, 1, 2, ones(Float64, 1,1,1)) isa NQCCalculators.RingPolymer_ClassicalFrictionModel_Cache
end

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
        @test NQCCalculators.get_derivative(cache, r)[1] !== NQCModels.derivative(model, r)
        @test NQCCalculators.get_eigen(cache, r).w !== eigvals(NQCModels.potential(model, r))
        @test abs.(NQCCalculators.get_eigen(cache, r).Z) !== abs.(eigvecs(NQCModels.potential(model, r)))
        @test NQCCalculators.get_adiabatic_derivative(cache, r)[1] !== cache.eigen.Z' * NQCModels.derivative(model, r)[1] * cache.eigen.Z
    end

    @testset "Dependent evaluation" begin
        #Change position value
        r = hcat(0.2)

        #Check that get_:quantities retrieves the wrong values
        @test NQCCalculators.get_potential(cache, r) !== NQCModels.potential(model, r)
        @test NQCCalculators.get_derivative(cache, r)[1] !== NQCModels.derivative(model, r)
        @test NQCCalculators.get_eigen(cache, r).w !== eigvals(NQCModels.potential(model, r))
        @test abs.(NQCCalculators.get_eigen(cache, r).Z) !== abs.(eigvecs(NQCModels.potential(model, r)))
        @test NQCCalculators.get_adiabatic_derivative(cache, r)[1] !== cache.eigen.Z' * NQCModels.derivative(model, r)[1] * cache.eigen.Z

		@debug "Cache info" cache
        #Update entries in the Cache to the new values
        NQCCalculators.update_cache!(cache, r)
		@debug "Cache info" cache

        # Check all the results
        @test cache.potential ≈ NQCModels.potential(model, r)
        @test cache.derivative[1] ≈ NQCModels.derivative(model, r)[1]
        @test cache.eigen.w ≈ eigvals(NQCModels.potential(model, r))
        @test abs.(cache.eigen.Z) ≈ abs.(eigvecs(NQCModels.potential(model, r)))
        @test isapprox(
            cache.adiabatic_derivative[1],
            cache.eigen.Z' * NQCModels.derivative(model, r)[1] * cache.eigen.Z
        )
    end

    @testset "Zero allocations" begin
        @test @allocated(NQCCalculators.get_potential(cache, r)) == 0
        @test @allocated(NQCCalculators.get_derivative(cache, r)) == 0
        @test @allocated(NQCCalculators.get_eigen(cache, r)) == 0
        @test @allocated(NQCCalculators.get_adiabatic_derivative(cache, r)) == 0
        @test @allocated(NQCCalculators.get_nonadiabatic_coupling(cache, r)) == 0
    end

    @testset "Explict Bath Model" begin
        model = NQCModels.AndersonHolstein(NQCModels.ErpenbeckThoss(;Γ=1.0), TrapezoidalRule(40, -1.0, 1.0))
        r = rand(1,1)
        cache = NQCCalculators.QuantumModel_Cache(model, 1, r)

        NQCCalculators.update_cache!(cache, r)

        @test cache.potential ≈ NQCModels.potential(model, r)
        @test cache.derivative ≈ NQCModels.derivative(model, r)
        @test cache.eigen.w ≈ eigvals(NQCModels.potential(model, r))
        @test abs.(cache.eigen.Z) ≈ abs.(eigvecs(NQCModels.potential(model, r)))
        @test isapprox(
            cache.adiabatic_derivative[1],
            cache.eigen.Z' * NQCModels.derivative(model, r)[1] * cache.eigen.Z
        )
    end
end

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
            @test NQCCalculators.get_eigen(cache, r)[i].w !== eigvals(NQCModels.potential(model, r[:,:,i]))
            @test abs.(NQCCalculators.get_eigen(cache, r)[i].Z) !== abs.(eigvecs(NQCModels.potential(model, r[:,:,i])))
            @test NQCCalculators.get_adiabatic_derivative(cache, r)[i] !== cache.eigen[i].Z' * NQCModels.derivative(model, r[:,:,i])[1] * cache.eigen[i].Z
        end
		@debug "Cache info" cache
        NQCCalculators.update_cache!(cache, r)
		@debug "Cache info" cache
        # Check that everything has been updated
        for i in NQCCalculators.beads(cache)
            @test NQCCalculators.get_derivative(cache, r)[i] == NQCModels.derivative(model, r[:,:,i])[1]
            @test NQCCalculators.get_eigen(cache, r)[i].w == eigvals(NQCModels.potential(model, r[:,:,i]))
            @test abs.(NQCCalculators.get_eigen(cache, r)[i].Z) == abs.(eigvecs(NQCModels.potential(model, r[:,:,i])))
            @test NQCCalculators.get_adiabatic_derivative(cache, r)[i] == cache.eigen[i].Z' * NQCModels.derivative(model, r[:,:,i])[1] * cache.eigen[i].Z
        end
    end

    @testset "Dependent evaluation" begin
        cache = reset_cache()

        NQCCalculators.update_cache!(cache, r)

        # Check all the results
        @test cache.potential ≈ [NQCModels.potential(model, r[:,:,i]) for i=1:10]
        @test cache.derivative ≈ [NQCModels.derivative(model, hcat(r[k,j,i]))[1] for k=1:1, j=1:1, i=1:10]
        @test [cache.eigen[i].w for i=1:10] ≈ eigvals.(cache.potential)
        @test [abs.(cache.eigen[i].Z) for i=1:10] ≈ [abs.(eigvecs(cache.potential[i])) for i=1:10]
        @test isapprox(
            cache.adiabatic_derivative[1],
            cache.eigen[1].Z' * NQCModels.derivative(model, r[:,:,1])[1] * cache.eigen[1].Z
        )
    end

    @testset "Centroid dependent evaluation" begin
        cache = reset_cache()
        r_centroid = RingPolymerArrays.get_centroid(r)


        NQCCalculators.update_cache!(cache, r)

        # Check all the results
        @test cache.centroid_potential ≈ NQCModels.potential(model, r_centroid)
        @test cache.centroid_derivative[1] ≈ NQCModels.derivative(model, r_centroid)[1]
        @test cache.centroid_eigen.w ≈ eigvals(NQCModels.potential(model, r_centroid))
        @test abs.(cache.centroid_eigen.Z) ≈ abs.(eigvecs(NQCModels.potential(model, r_centroid)))
        @test isapprox(
            cache.centroid_adiabatic_derivative[1],
            cache.centroid_eigen.Z' * NQCModels.derivative(model, r_centroid)[1] * cache.centroid_eigen.Z
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

@testset "RingPolymer_QuantumFrictionModel_Cache" begin
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
