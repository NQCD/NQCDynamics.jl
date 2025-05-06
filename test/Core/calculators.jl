using Test
using NQCDynamics
using NQCCalculators
using NQCModels
using LinearAlgebra: tr, Diagonal, eigvecs, eigvals
using RingPolymerArrays: RingPolymerArrays

# Allocation-based tests aren't great - change these values if there's a good reason to. 
# Eigenvalue allocations limit: Regardless of maths backend and Julia version, eigenvalues should take fewer allocations than this or something might be wrong
const MAX_EIGEN_ALLOCATIONS = 60000

@testset "General constructors" begin #reasonable basic tests to ensure that multiple dispatch is correctly working to construct the correct types of Calculators
    model = NQCModels.DoubleWell()
    @test NQCCalculators.Create_Cache(model, 1, Float64) isa NQCCalculators.SmallQuantumModel_Cache
    @test NQCCalculators.Create_Cache(model, 1, 2, Float64) isa NQCCalculators.RingPolymer_SmallQuantumModel_Cache

    model = NQCModels.Free()
    @test NQCCalculators.Create_Cache(model, 1, Float64) isa NQCCalculators.ClassicalModel_Cache
    @test NQCCalculators.Create_Cache(model, 1, 2, Float64) isa NQCCalculators.RingPolymer_ClassicalModel_Cache

    model = NQCModels.AndersonHolstein(MiaoSubotnik(;Γ=1.0), TrapezoidalRule(10, -1.0, 1.0))
    @test NQCCalculators.Create_Cache(model, 1, Float64) isa NQCCalculators.LargeQuantumModel_Cache

    model = NQCModels.CompositeFrictionModel(Free(), RandomFriction(1))
    @test NQCCalculators.Create_Cache(model, 1, Float64) isa NQCCalculators.Friction_Cache
    @test NQCCalculators.Create_Cache(model, 1, 2, Float64) isa NQCCalculators.RingPolymer_Friction_Cache
end

@testset "ClassicalModel_Cache" begin
    model = NQCModels.Free(3)
    cache = NQCCalculators.ClassicalModel_Cache{Float64}(model, 5) 
    r = rand(3, 5)

    @test @allocated(NQCCalculators.evaluate_potential!(cache, r)) == 0
    @test @allocated(NQCCalculators.evaluate_derivative!(cache, r)) == 0

    @test @allocated(NQCCalculators.get_potential(cache, r)) == 0
    @test @allocated(NQCCalculators.get_derivative(cache, r)) == 0

    #model = NQCModels.Morse(Float64)
    #cache = NQCCalculators.ClassicalModel_Cache{Float64}(model, 3)
end

@testset "RingPolymer_ClassicalModel_Cache" begin
    model = NQCModels.Free(3)
    cache = NQCCalculators.RingPolymer_ClassicalModel_Cache{Float64}(model, 5, 10) 
    r = rand(3, 5, 10)

    @test @allocated(NQCCalculators.evaluate_potential!(cache, r)) == 0
    @test @allocated(NQCCalculators.evaluate_derivative!(cache, r)) == 0

    @test @allocated(NQCCalculators.get_potential(cache, r)) == 0
    @test @allocated(NQCCalculators.get_derivative(cache, r)) == 0
end

@testset "SmallQuantumModel_Cache" begin
    model = NQCModels.DoubleWell()
    cache = NQCCalculators.SmallQuantumModel_Cache{Float64}(model, 1) 
    r = rand(1,1)
    quantities = (:potential, :derivative, :eigen, :adiabatic_derivative, :nonadiabatic_coupling)

    @testset "Potential evaluation" begin
        #check that the position is only changed in cache after get_potential() is called and not when potential() is called
        r = [0.1;;]
        true_potential = potential(model, r)

        @test cache.potential != true_potential
        @test getfield(cache, :potential).position != r

        NQCCalculators.get_potential(cache, r)

        @test cache.potential == true_potential

        # Check that the position for potential is the only one that has been updated. 
        for name in quantities
            if name !== :potential
                @test getfield(cache, name).position != r
            else
                @test getfield(cache, name).position == r
            end
        end
    end

    @testset "Dependent evaluation" begin
        r = [0.0;;]

        # Check the position fields are different
        for name in quantities
            @test getfield(cache, name).position != r
        end

        NQCCalculators.get_nonadiabatic_coupling(cache, r)

        # Check the position fields are correct
        for name in quantities
            @test getfield(cache, name).position == r
        end

        # Check all the results
        @test cache.potential ≈ potential(model, r)
        @test cache.derivative ≈ derivative(model, r)
        @test cache.eigen.values ≈ eigvals(potential(model, r))
        @test abs.(cache.eigen.vectors) ≈ abs.(eigvecs(potential(model, r)))
        @test isapprox(
            cache.adiabatic_derivative[1],
            cache.eigen.vectors' * derivative(model, r)[1] * cache.eigen.vectors
        )
    end

    @testset "Zero allocations" begin
        @test @allocated(NQCCalculators.evaluate_potential!(cache, r)) == 0
        @test @allocated(NQCCalculators.evaluate_derivative!(cache, r)) == 0
        @test @allocated(NQCCalculators.evaluate_eigen!(cache, r)) == 0
        @test @allocated(NQCCalculators.evaluate_adiabatic_derivative!(cache, r)) == 0
        @test @allocated(NQCCalculators.evaluate_nonadiabatic_coupling!(cache, r)) == 0
    end
end

@testset "RingPolymer_SmallQuantumModel_Cache" begin
    model = NQCModels.DoubleWell()

    reset_cache() = NQCCalculators.RingPolymer_SmallQuantumModel_Cache{Float64}(model, 1, 10)

    r = rand(1,1,10)
    r_centroid = rand(1,1)
    standard_quantities = (:potential, :derivative, :eigen, :adiabatic_derivative, :nonadiabatic_coupling)
    centroid_quantities = (:centroid, :centroid_potential, :centroid_derivative,
        :centroid_eigen, :centroid_adiabatic_derivative, :centroid_nonadiabatic_coupling)

    @testset "Potential evaluation" begin
        cache = reset_cache()
        true_potential = [potential(model, r[:,:,i]) for i=1:10]

        @test cache.potential != true_potential
        @test getfield(cache, :potential).position != r

        NQCCalculators.get_potential(cache, r)

        @test cache.potential == true_potential

        # Check that the position for potential is the only one that has been updated. 
        for name in standard_quantities
            if name !== :potential
                @test getfield(cache, name).position != r
            else
                @test getfield(cache, name).position == r
            end
        end
    end

    @testset "Dependent evaluation" begin
        cache = reset_cache()

        # Reset the position fields
        for name in standard_quantities
            getfield(cache, name).position .= NaN
        end

        NQCCalculators.get_nonadiabatic_coupling(cache, r)

        # Check the position fields are correct
        for name in standard_quantities
            @test getfield(cache, name).position == r
        end

        # Check all the results
        @test cache.potential ≈ [potential(model, r[:,:,i]) for i=1:10]
        @test cache.derivative ≈ [derivative(model, r[k,j,i]) for k=1:1, j=1:1, i=1:10]
        @test [cache.eigen[i].values for i=1:10] ≈ eigvals.(cache.potential)
        @test [abs.(cache.eigen[i].vectors) for i=1:10] ≈ [abs.(eigvecs(cache.potential[i])) for i=1:10]
        @test isapprox(
            cache.adiabatic_derivative[1],
            cache.eigen[1].vectors' * derivative(model, r[:,:,1])[1] * cache.eigen[1].vectors
        )
    end


    @testset "Centroid dependent evaluation" begin
        cache = reset_cache()
        r_centroid = RingPolymerArrays.get_centroid(r)

        # Check the position fields are different
        for name in centroid_quantities
            @test getfield(cache, name).position != r
        end

        NQCCalculators.get_centroid_nonadiabatic_coupling(cache, r)

        # Check the position fields are correct
        for name in centroid_quantities
            @test getfield(cache, name).position == r
        end

        # Check all the results
        @test cache.centroid_potential ≈ potential(model, r_centroid)
        @test cache.centroid_derivative ≈ derivative(model, r_centroid)
        @test cache.centroid_eigen.values ≈ eigvals(potential(model, r_centroid))
        @test abs.(cache.centroid_eigen.vectors) ≈ abs.(eigvecs(potential(model, r_centroid)))
        @test isapprox(
            cache.centroid_adiabatic_derivative[1],
            cache.centroid_eigen.vectors' * derivative(model, r_centroid)[1] * cache.centroid_eigen.vectors
        )
    end

    @testset "Extra quantities" begin
        cache = reset_cache()

        NQCCalculators.get_traceless_potential(cache, r)
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
        NQCCalculators.evaluate_potential!(cache, r)
        NQCCalculators.evaluate_derivative!(cache, r)
        NQCCalculators.evaluate_eigen!(cache, r)
        NQCCalculators.evaluate_adiabatic_derivative!(cache, r)
        NQCCalculators.evaluate_nonadiabatic_coupling!(cache, r)

        NQCCalculators.evaluate_traceless_potential!(cache, r)
        NQCCalculators.evaluate_V̄!(cache, r)
        NQCCalculators.evaluate_traceless_derivative!(cache, r)
        NQCCalculators.evaluate_D̄!(cache, r)
        NQCCalculators.evaluate_traceless_adiabatic_derivative!(cache, r)

        NQCCalculators.evaluate_centroid!(cache, r)
        NQCCalculators.evaluate_centroid_potential!(cache, r)
        NQCCalculators.evaluate_centroid_derivative!(cache, r)
        NQCCalculators.evaluate_centroid_eigen!(cache, r)
        NQCCalculators.evaluate_centroid_adiabatic_derivative!(cache, r)
        NQCCalculators.evaluate_centroid_nonadiabatic_coupling!(cache, r)

        cache = reset_cache()

        # Check for no allocations
        @test @allocated(NQCCalculators.evaluate_potential!(cache, r)) == 0
        @test @allocated(NQCCalculators.evaluate_derivative!(cache, r)) == 0
        @test @allocated(NQCCalculators.evaluate_eigen!(cache, r)) == 0
        @test @allocated(NQCCalculators.evaluate_adiabatic_derivative!(cache, r)) == 0
        @test @allocated(NQCCalculators.evaluate_nonadiabatic_coupling!(cache, r)) == 0

        @test @allocated(NQCCalculators.evaluate_traceless_potential!(cache, r)) == 0
        @test @allocated(NQCCalculators.evaluate_V̄!(cache, r)) == 0
        @test @allocated(NQCCalculators.evaluate_traceless_derivative!(cache, r)) == 0
        @test @allocated(NQCCalculators.evaluate_D̄!(cache, r)) == 0
        @test @allocated(NQCCalculators.evaluate_traceless_adiabatic_derivative!(cache, r)) == 0

        @test @allocated(NQCCalculators.evaluate_centroid!(cache, r)) == 0
        @test @allocated(NQCCalculators.evaluate_centroid_potential!(cache, r)) == 48 # not sure why this isn't zero
        @test @allocated(NQCCalculators.evaluate_centroid_derivative!(cache, r)) == 0
        @test @allocated(NQCCalculators.evaluate_centroid_eigen!(cache, r)) == 0
        @test @allocated(NQCCalculators.evaluate_centroid_adiabatic_derivative!(cache, r)) == 0
        @test @allocated(NQCCalculators.evaluate_centroid_nonadiabatic_coupling!(cache, r)) == 0
    end
end

@testset "LargeQuantumModel_Cache" begin
    model = NQCModels.AndersonHolstein(MiaoSubotnik(;Γ=1.0), TrapezoidalRule(40, -1.0, 1.0))
    cache = NQCCalculators.LargeQuantumModel_Cache{Float64}(model, 1) 
    r = rand(1,1)

    NQCCalculators.get_nonadiabatic_coupling(cache, r)
    @test all(x->x==1, values(cache.stats))

    NQCCalculators.evaluate_potential!(cache, r)
    NQCCalculators.evaluate_derivative!(cache, r)
    NQCCalculators.evaluate_eigen!(cache, r)
    NQCCalculators.evaluate_adiabatic_derivative!(cache, r)
    NQCCalculators.evaluate_nonadiabatic_coupling!(cache, r)

    @test @allocated(NQCCalculators.evaluate_potential!(cache, r)) == 0
    @test @allocated(NQCCalculators.evaluate_derivative!(cache, r)) == 0
    eigen_allocations = @allocated(NQCCalculators.evaluate_eigen!(cache, r))
    global allocs_LargeDiabatic = eigen_allocations
    @debug "LargeQuantumModel_Cache: $(eigen_allocations) needed for eigenvalues"
    @test eigen_allocations < MAX_EIGEN_ALLOCATIONS
    @test @allocated(NQCCalculators.evaluate_adiabatic_derivative!(cache, r)) == 0
    @test @allocated(NQCCalculators.evaluate_nonadiabatic_coupling!(cache, r)) == 0
end

@testset "RingPolymer_LargeQuantumModel_Cache" begin
    model = NQCModels.AndersonHolstein(MiaoSubotnik(;Γ=1.0), TrapezoidalRule(40, -1.0, 1.0))
    cache = NQCCalculators.RingPolymer_LargeQuantumModel_Cache{Float64}(model, 1, 10) 
    r = rand(1, 1, 10)

    NQCCalculators.get_nonadiabatic_coupling(cache, r)
    NQCCalculators.get_centroid_nonadiabatic_coupling(cache, r)
    @test all(x->x==1, values(cache.stats))

    NQCCalculators.evaluate_potential!(cache, r)
    NQCCalculators.evaluate_derivative!(cache, r)
    NQCCalculators.evaluate_eigen!(cache, r)
    NQCCalculators.evaluate_adiabatic_derivative!(cache, r)
    NQCCalculators.evaluate_nonadiabatic_coupling!(cache, r)

    @test @allocated(NQCCalculators.evaluate_potential!(cache, r)) == 0
    @test @allocated(NQCCalculators.evaluate_derivative!(cache, r)) == 0
    @test @allocated(NQCCalculators.evaluate_eigen!(cache, r)) == 10 * allocs_LargeDiabatic
    @test @allocated(NQCCalculators.evaluate_adiabatic_derivative!(cache, r)) == 0
    @test @allocated(NQCCalculators.evaluate_nonadiabatic_coupling!(cache, r)) == 0
end

@testset "Friction_Cache" begin
    model = CompositeFrictionModel(Free(), RandomFriction(1))
    cache = NQCCalculators.Friction_Cache{Float64}(model, 1) 
    r = rand(1,1)

    NQCCalculators.get_potential(cache, r)
    NQCCalculators.get_derivative(cache, r)
    NQCCalculators.get_friction(cache, r)
    @test all(x->x==1, values(cache.stats))

    NQCCalculators.evaluate_potential!(cache, r)
    NQCCalculators.evaluate_derivative!(cache, r)
    NQCCalculators.evaluate_friction!(cache, r)
    @test @allocated(NQCCalculators.evaluate_potential!(cache, r)) == 0
    @test @allocated(NQCCalculators.evaluate_derivative!(cache, r)) == 0
    global allocs_Friction = @allocated(NQCCalculators.evaluate_friction!(cache, r)) 
    @debug "Friction_Cache.evaluate_friction made $(allocs_Friction) allocations. This should be around 250."
    @test allocs_Friction ≤ 250
end

@testset "RingPolymer_Friction_Cache" begin
    model = CompositeFrictionModel(Free(), RandomFriction(1))
    cache = NQCCalculators.RingPolymer_Friction_Cache{Float64}(model, 1, 10)
    r = rand(1,1,10)

    NQCCalculators.get_potential(cache, r)
    NQCCalculators.get_derivative(cache, r)
    NQCCalculators.get_friction(cache, r)
    @test all(x->x==1, values(cache.stats))

    NQCCalculators.evaluate_potential!(cache, r)
    NQCCalculators.evaluate_derivative!(cache, r)
    NQCCalculators.evaluate_friction!(cache, r)
    @test @allocated(NQCCalculators.evaluate_potential!(cache, r)) == 0
    @test @allocated(NQCCalculators.evaluate_derivative!(cache, r)) == 0
    @test @allocated(NQCCalculators.evaluate_friction!(cache, r)) == 10 * allocs_Friction
end
