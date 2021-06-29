using Test
using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.Calculators

@testset "Adiabatic" begin
    model = NonadiabaticModels.Free()
    @test Calculators.AdiabaticCalculator{Float64}(model, 3, 5) isa Calculators.AdiabaticCalculator
    @test Calculators.RingPolymerAdiabaticCalculator{Float64}(model, 3, 5, 10) isa Calculators.RingPolymerAdiabaticCalculator

    calc = Calculators.AdiabaticCalculator{Float64}(model, 3, 5) 
    Calculators.evaluate_potential!(calc, rand(3, 5))
    Calculators.evaluate_derivative!(calc, rand(3, 5))

    calc = Calculators.RingPolymerAdiabaticCalculator{Float64}(model, 3, 5, 10) 
    Calculators.evaluate_potential!(calc, rand(3, 5, 10))
    Calculators.evaluate_derivative!(calc, rand(3, 5, 10))
end

@testset "Diabatic" begin
    model = NonadiabaticModels.DoubleWell()
    @test Calculators.DiabaticCalculator{Float64}(model, 3, 5) isa Calculators.DiabaticCalculator
    @test Calculators.RingPolymerDiabaticCalculator{Float64}(model, 3, 5, 10) isa Calculators.RingPolymerDiabaticCalculator

    calc = Calculators.DiabaticCalculator{Float64}(model, 1, 1) 
    Calculators.evaluate_potential!(calc, rand(1, 1))
    Calculators.evaluate_derivative!(calc, rand(1, 1))
    Calculators.eigen!(calc)
    Calculators.transform_derivative!(calc)
    Calculators.evaluate_nonadiabatic_coupling!(calc)
    Calculators.update_electronics!(calc, rand(1, 1))
    
    calc = Calculators.RingPolymerDiabaticCalculator{Float64}(model, 1, 1, 10)
    Calculators.evaluate_potential!(calc, rand(1, 1, 10))
    Calculators.evaluate_centroid_potential!(calc, rand(1, 1, 10))
    Calculators.evaluate_derivative!(calc, rand(1, 1, 10))
    Calculators.evaluate_centroid_derivative!(calc, rand(1, 1, 10))
    Calculators.eigen!(calc)
    Calculators.transform_derivative!(calc)
    Calculators.evaluate_nonadiabatic_coupling!(calc)
    Calculators.update_electronics!(calc, rand(1, 1, 10))
end

@testset "LargeDiabatic" begin
    model = NonadiabaticModels.MiaoSubotnik()
    @test Calculators.LargeDiabaticCalculator{Float64}(model, 3, 5) isa Calculators.LargeDiabaticCalculator

    calc = Calculators.LargeDiabaticCalculator{Float64}(model, 1, 1) 
    Calculators.evaluate_potential!(calc, rand(1, 1))
    Calculators.evaluate_derivative!(calc, rand(1, 1))
    Calculators.eigen!(calc)
    Calculators.transform_derivative!(calc)
    Calculators.evaluate_nonadiabatic_coupling!(calc)
    Calculators.update_electronics!(calc, rand(1, 1))
end

@testset "General constructors" begin
    model = NonadiabaticModels.DoubleWell()
    @test Calculators.Calculator(model, 1, 1) isa Calculators.DiabaticCalculator
    @test Calculators.Calculator(model, 1, 1, 2) isa Calculators.RingPolymerDiabaticCalculator
    model = NonadiabaticModels.Free()
    @test Calculators.Calculator(model, 1, 1) isa Calculators.AdiabaticCalculator
    @test Calculators.Calculator(model, 1, 1, 2) isa Calculators.RingPolymerAdiabaticCalculator
    model = NonadiabaticModels.MiaoSubotnik()
    @test Calculators.Calculator(model, 1, 1) isa Calculators.LargeDiabaticCalculator
end
