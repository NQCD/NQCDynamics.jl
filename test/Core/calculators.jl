using Test
using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.Calculators
using LinearAlgebra: tr, Diagonal

@testset "Adiabatic" begin
    model = NonadiabaticModels.Free(3)
    @test Calculators.AdiabaticCalculator{Float64}(model, 5) isa Calculators.AdiabaticCalculator
    @test Calculators.RingPolymerAdiabaticCalculator{Float64}(model, 5, 10) isa Calculators.RingPolymerAdiabaticCalculator

    calc = Calculators.AdiabaticCalculator{Float64}(model, 5) 
    Calculators.evaluate_potential!(calc, rand(3, 5))
    Calculators.evaluate_derivative!(calc, rand(3, 5))

    calc = Calculators.RingPolymerAdiabaticCalculator{Float64}(model, 5, 10) 
    Calculators.evaluate_potential!(calc, rand(3, 5, 10))
    Calculators.evaluate_derivative!(calc, rand(3, 5, 10))
end

@testset "Diabatic" begin
    model = NonadiabaticModels.DoubleWell()
    @test Calculators.DiabaticCalculator{Float64}(model, 5) isa Calculators.DiabaticCalculator
    @test Calculators.RingPolymerDiabaticCalculator{Float64}(model, 5, 10) isa Calculators.RingPolymerDiabaticCalculator

    calc = Calculators.DiabaticCalculator{Float64}(model, 1) 
    Calculators.evaluate_potential!(calc, rand(1, 1))
    Calculators.evaluate_derivative!(calc, rand(1, 1))
    Calculators.eigen!(calc)
    Calculators.transform_derivative!(calc)
    Calculators.evaluate_nonadiabatic_coupling!(calc)
    Calculators.update_electronics!(calc, rand(1, 1))
    
    calc = Calculators.RingPolymerDiabaticCalculator{Float64}(model, 1, 10)
    Calculators.evaluate_potential!(calc, rand(1, 1, 10))
    Calculators.evaluate_centroid_potential!(calc, rand(1, 1))

    Calculators.evaluate_V̄!(calc)
    @test calc.V̄[1] ≈ tr(calc.potential[1]) / nstates(model)
    @test (@allocated Calculators.evaluate_V̄!(calc)) == 0

    Calculators.evaluate_traceless_potential!(calc)
    @test calc.traceless_potential[1] ≈ calc.potential[1] .- Diagonal(fill(calc.V̄[1], nstates(model)))
    @test tr(calc.traceless_potential[1]) ≈ 0 atol=eps(Float64)
    @test (@allocated Calculators.evaluate_traceless_potential!(calc)) == 0

    Calculators.evaluate_derivative!(calc, rand(1, 1, 10))
    Calculators.evaluate_centroid_derivative!(calc, rand(1, 1))

    Calculators.evaluate_D̄!(calc)
    @test calc.D̄[1] ≈ tr(calc.derivative[1]) / nstates(model)
    @test (@allocated Calculators.evaluate_D̄!(calc)) == 0

    Calculators.evaluate_traceless_derivative!(calc)
    @test calc.traceless_derivative[1] ≈ calc.derivative[1] .- Diagonal(fill(calc.D̄[1], nstates(model)))
    @test tr(calc.traceless_derivative[1]) ≈ 0 atol=eps(Float64)
    @test (@allocated Calculators.evaluate_traceless_derivative!(calc)) == 0

    Calculators.evaluate_traceless_adiabatic_derivative!(calc)
    @test calc.traceless_adiabatic_derivative[1] ≈ calc.adiabatic_derivative[1] .- Diagonal(fill(calc.D̄[1], nstates(model)))
    @test (@allocated Calculators.evaluate_traceless_adiabatic_derivative!(calc)) == 0

    Calculators.eigen!(calc)
    Calculators.transform_derivative!(calc)
    Calculators.evaluate_nonadiabatic_coupling!(calc)
    Calculators.update_electronics!(calc, rand(1, 1, 10))
end

@testset "LargeDiabatic" begin
    model = NonadiabaticModels.MiaoSubotnik()
    @test Calculators.LargeDiabaticCalculator{Float64}(model, 5) isa Calculators.LargeDiabaticCalculator

    calc = Calculators.LargeDiabaticCalculator{Float64}(model, 1) 
    Calculators.evaluate_potential!(calc, rand(1, 1))
    Calculators.evaluate_derivative!(calc, rand(1, 1))
    Calculators.eigen!(calc)
    Calculators.transform_derivative!(calc)
    Calculators.evaluate_nonadiabatic_coupling!(calc)
    Calculators.update_electronics!(calc, rand(1, 1))
end

@testset "General constructors" begin
    model = NonadiabaticModels.DoubleWell()
    @test Calculators.Calculator(model, 1, Float64) isa Calculators.DiabaticCalculator
    @test Calculators.Calculator(model, 1, 2, Float64) isa Calculators.RingPolymerDiabaticCalculator
    model = NonadiabaticModels.Free()
    @test Calculators.Calculator(model, 1, Float64) isa Calculators.AdiabaticCalculator
    @test Calculators.Calculator(model, 1, 2, Float64) isa Calculators.RingPolymerAdiabaticCalculator
    model = NonadiabaticModels.MiaoSubotnik()
    @test Calculators.Calculator(model, 1, Float64) isa Calculators.LargeDiabaticCalculator
end

# for field ∈ (:potential, :derivative)
#     @testset "get_$field" begin
#         get_field = Symbol(:get, :_, field)
#         evaluate_field = Symbol(:evaluate, :_, field, :!)
#         eval(quote
#             model = NonadiabaticModels.DoubleWell()
#             calc = Calculators.DiabaticCalculator{Float64}(model, 1)
#             r = hcat(1.0)
#             @test !(Calculators.position(calc) ≈ r)
#             Calculators.$get_field(calc, r)
#             @test Calculators.position(calc) ≈ r
#             val1 = Calculators.$get_field(calc, r) 
#             val2 = Calculators.$evaluate_field(calc, r)
#             @test val1 === val2 === calc.$field
#         end)
#     end
# end