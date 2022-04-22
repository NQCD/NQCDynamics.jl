using Test
using NQCDynamics
using NQCDynamics.Calculators
using LinearAlgebra: tr, Diagonal, eigvecs, eigvals
using RingPolymerArrays: RingPolymerArrays

@testset "General constructors" begin
    model = NQCModels.DoubleWell()
    @test Calculators.Calculator(model, 1, Float64) isa Calculators.DiabaticCalculator
    @test Calculators.Calculator(model, 1, 2, Float64) isa Calculators.RingPolymerDiabaticCalculator
    model = NQCModels.Free()
    @test Calculators.Calculator(model, 1, Float64) isa Calculators.AdiabaticCalculator
    @test Calculators.Calculator(model, 1, 2, Float64) isa Calculators.RingPolymerAdiabaticCalculator
    model = NQCModels.MiaoSubotnik()
    @test Calculators.Calculator(model, 1, Float64) isa Calculators.LargeDiabaticCalculator
end

@testset "AdiabaticCalculator" begin
    model = NQCModels.Free(3)
    calc = Calculators.AdiabaticCalculator{Float64}(model, 5) 
    r = rand(3, 5)

    @test @allocated(Calculators.evaluate_potential!(calc, r)) == 0
    @test @allocated(Calculators.evaluate_derivative!(calc, r)) == 0

    @test @allocated(Calculators.get_potential(calc, r)) == 0
    @test @allocated(Calculators.get_derivative(calc, r)) == 0
end

@testset "RingPolymerAdiabaticCalculator" begin
    model = NQCModels.Free(3)
    calc = Calculators.RingPolymerAdiabaticCalculator{Float64}(model, 5, 10) 
    r = rand(3, 5, 10)

    @test @allocated(Calculators.evaluate_potential!(calc, r)) == 0
    @test @allocated(Calculators.evaluate_derivative!(calc, r)) == 0

    @test @allocated(Calculators.get_potential(calc, r)) == 0
    @test @allocated(Calculators.get_derivative(calc, r)) == 0
end

@testset "DiabaticCalculator" begin
    model = NQCModels.DoubleWell()
    calc = Calculators.DiabaticCalculator{Float64}(model, 1) 
    r = rand(1,1)
    quantities = (:potential, :derivative, :eigen, :adiabatic_derivative, :nonadiabatic_coupling)

    @testset "Potential evaluation" begin
        r = [0.1;;]
        true_potential = potential(model, r)

        @test calc.potential != true_potential
        @test getfield(calc, :potential).position != r

        Calculators.get_potential(calc, r)

        @test calc.potential == true_potential

        # Check that the position for potential is the only one that has been updated. 
        for name in quantities
            if name !== :potential
                @test getfield(calc, name).position != r
            else
                @test getfield(calc, name).position == r
            end
        end
    end

    @testset "Dependent evaluation" begin
        r = [0.0;;]

        # Check the position fields are different
        for name in quantities
            @test getfield(calc, name).position != r
        end

        Calculators.get_nonadiabatic_coupling(calc, r)

        # Check the position fields are correct
        for name in quantities
            @test getfield(calc, name).position == r
        end

        # Check all the results
        @test calc.potential ≈ potential(model, r)
        @test calc.derivative ≈ derivative(model, r)
        @test calc.eigen.values ≈ eigvals(potential(model, r))
        @test abs.(calc.eigen.vectors) ≈ abs.(eigvecs(potential(model, r)))
        @test isapprox(
            calc.adiabatic_derivative[1],
            calc.eigen.vectors' * derivative(model, r)[1] * calc.eigen.vectors
        )
    end

    @testset "Zero allocations" begin
        @test @allocated(Calculators.evaluate_potential!(calc, r)) == 0
        @test @allocated(Calculators.evaluate_derivative!(calc, r)) == 0
        @test @allocated(Calculators.evaluate_eigen!(calc, r)) == 0
        @test @allocated(Calculators.evaluate_adiabatic_derivative!(calc, r)) == 0
        @test @allocated(Calculators.evaluate_nonadiabatic_coupling!(calc, r)) == 0
    end
end

@testset "RingPolymerDiabaticCalculator" begin
    model = NQCModels.DoubleWell()

    reset_calc() = Calculators.RingPolymerDiabaticCalculator{Float64}(model, 1, 10)

    r = rand(1,1,10)
    r_centroid = rand(1,1)
    standard_quantities = (:potential, :derivative, :eigen, :adiabatic_derivative, :nonadiabatic_coupling)
    centroid_quantities = (:centroid, :centroid_potential, :centroid_derivative,
        :centroid_eigen, :centroid_adiabatic_derivative, :centroid_nonadiabatic_coupling)

    @testset "Potential evaluation" begin
        calc = reset_calc()
        true_potential = [potential(model, r[:,:,i]) for i=1:10]

        @test calc.potential != true_potential
        @test getfield(calc, :potential).position != r

        Calculators.get_potential(calc, r)

        @test calc.potential == true_potential

        # Check that the position for potential is the only one that has been updated. 
        for name in standard_quantities
            if name !== :potential
                @test getfield(calc, name).position != r
            else
                @test getfield(calc, name).position == r
            end
        end
    end

    @testset "Dependent evaluation" begin
        calc = reset_calc()

        # Reset the position fields
        for name in standard_quantities
            getfield(calc, name).position .= NaN
        end

        Calculators.get_nonadiabatic_coupling(calc, r)

        # Check the position fields are correct
        for name in standard_quantities
            @test getfield(calc, name).position == r
        end

        # Check all the results
        @test calc.potential ≈ [potential(model, r[:,:,i]) for i=1:10]
        @test calc.derivative ≈ [derivative(model, r[k,j,i]) for k=1:1, j=1:1, i=1:10]
        @test [calc.eigen[i].values for i=1:10] ≈ eigvals.(calc.potential)
        @test [abs.(calc.eigen[i].vectors) for i=1:10] ≈ [abs.(eigvecs(calc.potential[i])) for i=1:10]
        @test isapprox(
            calc.adiabatic_derivative[1],
            calc.eigen[1].vectors' * derivative(model, r[:,:,1])[1] * calc.eigen[1].vectors
        )
    end


    @testset "Centroid dependent evaluation" begin
        calc = reset_calc()
        r_centroid = RingPolymerArrays.get_centroid(r)

        # Check the position fields are different
        for name in centroid_quantities
            @test getfield(calc, name).position != r
        end

        Calculators.get_centroid_nonadiabatic_coupling(calc, r)

        # Check the position fields are correct
        for name in centroid_quantities
            @test getfield(calc, name).position == r
        end

        # Check all the results
        @test calc.centroid_potential ≈ potential(model, r_centroid)
        @test calc.centroid_derivative ≈ derivative(model, r_centroid)
        @test calc.centroid_eigen.values ≈ eigvals(potential(model, r_centroid))
        @test abs.(calc.centroid_eigen.vectors) ≈ abs.(eigvecs(potential(model, r_centroid)))
        @test isapprox(
            calc.centroid_adiabatic_derivative[1],
            calc.centroid_eigen.vectors' * derivative(model, r_centroid)[1] * calc.centroid_eigen.vectors
        )
    end

    @testset "Extra quantities" begin
        calc = reset_calc()

        Calculators.get_traceless_potential(calc, r)
        Calculators.get_traceless_derivative(calc, r)
        Calculators.get_traceless_adiabatic_derivative(calc, r)

        @test calc.V̄[1] ≈ tr(calc.potential[1]) / nstates(model)
        @test tr(calc.traceless_potential[1]) ≈ 0 atol=eps(Float64)
        @test calc.D̄[1] ≈ tr(calc.derivative[1]) / nstates(model)
        @test tr(calc.traceless_derivative[1]) ≈ 0 atol=eps(Float64)
        @test calc.traceless_adiabatic_derivative[1] ≈ calc.adiabatic_derivative[1] .- Diagonal(fill(calc.D̄[1], nstates(model)))
    end

    @testset "Zero allocations" begin

        calc = reset_calc()

        # Call all functions to ensure they have been compiled
        Calculators.evaluate_potential!(calc, r)
        Calculators.evaluate_derivative!(calc, r)
        Calculators.evaluate_eigen!(calc, r)
        Calculators.evaluate_adiabatic_derivative!(calc, r)
        Calculators.evaluate_nonadiabatic_coupling!(calc, r)

        Calculators.evaluate_traceless_potential!(calc, r)
        Calculators.evaluate_V̄!(calc, r)
        Calculators.evaluate_traceless_derivative!(calc, r)
        Calculators.evaluate_D̄!(calc, r)
        Calculators.evaluate_traceless_adiabatic_derivative!(calc, r)

        Calculators.evaluate_centroid!(calc, r)
        Calculators.evaluate_centroid_potential!(calc, r)
        Calculators.evaluate_centroid_derivative!(calc, r)
        Calculators.evaluate_centroid_eigen!(calc, r)
        Calculators.evaluate_centroid_adiabatic_derivative!(calc, r)
        Calculators.evaluate_centroid_nonadiabatic_coupling!(calc, r)

        calc = reset_calc()

        # Check for no allocations
        @test @allocated(Calculators.evaluate_potential!(calc, r)) == 0
        @test @allocated(Calculators.evaluate_derivative!(calc, r)) == 0
        @test @allocated(Calculators.evaluate_eigen!(calc, r)) == 0
        @test @allocated(Calculators.evaluate_adiabatic_derivative!(calc, r)) == 0
        @test @allocated(Calculators.evaluate_nonadiabatic_coupling!(calc, r)) == 0

        @test @allocated(Calculators.evaluate_traceless_potential!(calc, r)) == 0
        @test @allocated(Calculators.evaluate_V̄!(calc, r)) == 0
        @test @allocated(Calculators.evaluate_traceless_derivative!(calc, r)) == 0
        @test @allocated(Calculators.evaluate_D̄!(calc, r)) == 0
        @test @allocated(Calculators.evaluate_traceless_adiabatic_derivative!(calc, r)) == 0

        @test @allocated(Calculators.evaluate_centroid!(calc, r)) == 0
        @test @allocated(Calculators.evaluate_centroid_potential!(calc, r)) == 48 # not sure why this isn't zero
        @test @allocated(Calculators.evaluate_centroid_derivative!(calc, r)) == 0
        @test @allocated(Calculators.evaluate_centroid_eigen!(calc, r)) == 0
        @test @allocated(Calculators.evaluate_centroid_adiabatic_derivative!(calc, r)) == 0
        @test @allocated(Calculators.evaluate_centroid_nonadiabatic_coupling!(calc, r)) == 0
    end
end

@testset "LargeDiabaticCalculator" begin
    model = NQCModels.MiaoSubotnik()
    calc = Calculators.LargeDiabaticCalculator{Float64}(model, 1) 
    r = rand(1,1)

    Calculators.get_nonadiabatic_coupling(calc, r)
    @test all(x->x==1, values(calc.stats))

    Calculators.evaluate_potential!(calc, r)
    Calculators.evaluate_derivative!(calc, r)
    Calculators.evaluate_eigen!(calc, r)
    Calculators.evaluate_adiabatic_derivative!(calc, r)
    Calculators.evaluate_nonadiabatic_coupling!(calc, r)

    @test @allocated(Calculators.evaluate_potential!(calc, r)) == 0
    @test @allocated(Calculators.evaluate_derivative!(calc, r)) == 0
    @test @allocated(Calculators.evaluate_eigen!(calc, r)) == 56896 # nonzero due to eigenroutines
    @test @allocated(Calculators.evaluate_adiabatic_derivative!(calc, r)) == 0
    @test @allocated(Calculators.evaluate_nonadiabatic_coupling!(calc, r)) == 0
end

@testset "RingPolymerLargeDiabaticCalculator" begin
    model = NQCModels.MiaoSubotnik()
    calc = Calculators.RingPolymerLargeDiabaticCalculator{Float64}(model, 1, 10) 
    r = rand(1, 1, 10)

    Calculators.get_nonadiabatic_coupling(calc, r)
    @test all(x->x==1, values(calc.stats))

    Calculators.evaluate_potential!(calc, r)
    Calculators.evaluate_derivative!(calc, r)
    Calculators.evaluate_eigen!(calc, r)
    Calculators.evaluate_adiabatic_derivative!(calc, r)
    Calculators.evaluate_nonadiabatic_coupling!(calc, r)

    @test @allocated(Calculators.evaluate_potential!(calc, r)) == 0
    @test @allocated(Calculators.evaluate_derivative!(calc, r)) == 0
    @test @allocated(Calculators.evaluate_eigen!(calc, r)) == 568960 # nonzero due to eigenroutines
    @test @allocated(Calculators.evaluate_adiabatic_derivative!(calc, r)) == 0
    @test @allocated(Calculators.evaluate_nonadiabatic_coupling!(calc, r)) == 0
end

@testset "FrictionCalculator" begin
    model = CompositeFrictionModel(Free(), RandomFriction(1))
    calc = Calculators.FrictionCalculator{Float64}(model, 1) 
    r = rand(1,1)

    Calculators.get_potential(calc, r)
    Calculators.get_derivative(calc, r)
    Calculators.get_friction(calc, r)
    @test all(x->x==1, values(calc.stats))

    Calculators.evaluate_potential!(calc, r)
    Calculators.evaluate_derivative!(calc, r)
    Calculators.evaluate_friction!(calc, r)
    @test @allocated(Calculators.evaluate_potential!(calc, r)) == 0
    @test @allocated(Calculators.evaluate_derivative!(calc, r)) == 0
    @test @allocated(Calculators.evaluate_friction!(calc, r)) == 192
end

@testset "RingPolymerFrictionCalculator" begin
    model = CompositeFrictionModel(Free(), RandomFriction(1))
    calc = Calculators.RingPolymerFrictionCalculator{Float64}(model, 1, 10)
    r = rand(1,1,10)

    Calculators.get_potential(calc, r)
    Calculators.get_derivative(calc, r)
    Calculators.get_friction(calc, r)
    @test all(x->x==1, values(calc.stats))

    Calculators.evaluate_potential!(calc, r)
    Calculators.evaluate_derivative!(calc, r)
    Calculators.evaluate_friction!(calc, r)
    @test @allocated(Calculators.evaluate_potential!(calc, r)) == 0
    @test @allocated(Calculators.evaluate_derivative!(calc, r)) == 0
    @test @allocated(Calculators.evaluate_friction!(calc, r)) == 1920
end
