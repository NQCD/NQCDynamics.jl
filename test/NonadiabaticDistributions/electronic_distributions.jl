using Test
using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics: Calculators

@testset "fill_density!" begin
    @testset "SingleState" for state ∈ (1, 2)
        electronics = SingleState(state)
        density = zeros(2,2)
        NonadiabaticDistributions.fill_density!(density, electronics)
        @test density[state, state] == 1
    end

    @testset "ElectronicPopulation" begin
        population = [0.5, 0.5]
        electronics = ElectronicPopulation(population)
        density = zeros(2, 2)
        NonadiabaticDistributions.fill_density!(density, electronics)
        @test density[1, 1] == density[2, 2] == 0.5
    end
end

@testset "initialise_density_matrix" begin
    calc = Calculators.Calculator(DoubleWell(), 1, Float64)
    electronics = SingleState(1)
    density = NonadiabaticDistributions.initialise_density_matrix(electronics, calc)
    @test density == [1 0; 0 0]
end

rs = (fill(0.0, 1, 1), fill(0.0, 1, 1, 1))
calculators = (
    Calculators.Calculator(DoubleWell(), 1, Float64),
    Calculators.Calculator(DoubleWell(), 1, 1, Float64),
    )
@testset "transform_density!" for (r, calc) ∈ zip(rs, calculators)
    density = [1.0 0; 0 0]
    NonadiabaticDistributions.transform_density!(density, calc, r, :to_adiabatic)
    @test density ≈ [0.5 0.5; 0.5 0.5]
    NonadiabaticDistributions.transform_density!(density, calc, r, :to_diabatic)
    @test density ≈ [1.0 0; 0 0]
    @test_throws ArgumentError NonadiabaticDistributions.transform_density!(density, calc, r, :blah)
end

@testset "initialise_adiabatic_density_matrix" for (r, calc) ∈ zip(rs, calculators)
    @testset "Adiabatic" begin
        electronics = SingleState(1, Adiabatic())
        density = NonadiabaticDistributions.initialise_adiabatic_density_matrix(electronics, calc, r)
        @test density[1,1] == 1
    end

    @testset "Diabatic" begin
        electronics = SingleState(1, Diabatic())
        density = NonadiabaticDistributions.initialise_adiabatic_density_matrix(electronics, calc, r)
        @test all(density .≈ 0.5)
    end
end
