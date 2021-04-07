using Test
using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.InitialConditions.QuantisedDiatomic
using LinearAlgebra: norm

atoms = Atoms([:N, :O])
model = Models.DiatomicHarmonic()
sim = Simulation(atoms, model, Dynamics.Classical())

@testset "calculate_diatomic_energy" begin
    @test QuantisedDiatomic.calculate_diatomic_energy(sim, 1.0, normal_vector=rand(3)) ≈ 0.0 atol=1e-3
    @test QuantisedDiatomic.calculate_diatomic_energy(sim, 2.0, normal_vector=rand(3)) ≈ 0.5
end

@testset "calculate_force_constant" begin
    @test QuantisedDiatomic.calculate_force_constant(sim) ≈ 1
    model2 = Models.DiatomicHarmonic(r₀=6.3)
    sim2 = Simulation(atoms, model2, Dynamics.Classical())
    @test QuantisedDiatomic.calculate_force_constant(sim2, bond_length=4.0) ≈ 1
end

@testset "subtract_centre_of_mass!" begin
    atoms = Atoms([:O, :O])
    r = rand(3, 2)
    r = QuantisedDiatomic.subtract_centre_of_mass(r, atoms.masses)
    @test r[:,1] ≈ -r[:,2]
end

@testset "apply_random_rotation!" begin
    x = fill(5.0, 3, 2)
    y = fill(-4.0, 3, 2)
    @test all(x .== 5)
    @test all(y .== -4)
    QuantisedDiatomic.apply_random_rotation!(x, y)
    @test all(x .!= 5)
    @test all(y .!= -4)
    for col in eachcol(x)
        @test norm(col) ≈ norm([5, 5, 5])
    end
    for col in eachcol(y)
        @test norm(col) ≈ norm([-4, -4, -4])
    end
end


@testset "generate and check results" begin
    numbers = [rand(1:50), rand(1:50)]
    configs = generate_configurations(sim, numbers[1], numbers[2]; samples=10)
    for config in configs
        ν, J = quantise_diatomic(sim, config...)
        @test ν ≈ numbers[1] rtol=1e-1
        @test J ≈ numbers[2] rtol=1e-1
    end
end

