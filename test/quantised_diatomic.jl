using Test
using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.InitialConditions.QuantisedDiatomic
using LinearAlgebra: norm
using Unitful
using UnitfulAtomic

atoms = Atoms([:N, :O])
model = NonadiabaticModels.DiatomicHarmonic()
sim = Simulation(atoms, model, Dynamics.Classical())

@testset "calculate_diatomic_energy" begin
    @test QuantisedDiatomic.calculate_diatomic_energy(sim, 1.0, normal_vector=rand(3)) ≈ 0.0 atol=1e-3
    @test QuantisedDiatomic.calculate_diatomic_energy(sim, 2.0, normal_vector=rand(3)) ≈ 0.5
end

@testset "calculate_force_constant" begin
    @test QuantisedDiatomic.calculate_force_constant(sim) ≈ 1
    model2 = NonadiabaticModels.DiatomicHarmonic(r₀=6.3)
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
    numbers = [2, 5]
    configs = generate_configurations(sim, numbers[1], numbers[2]; samples=10, translational_energy=1u"eV")
    for config in configs
        ν, J = quantise_diatomic(sim, config...)
        @test ν ≈ numbers[1] rtol=1e-1
        @test J ≈ numbers[2] rtol=1e-1
    end
end

@testset "position_above_surface!" begin
    r = zeros(3, 2)
    height = 10
    direction = [0, 1, 0]
    QuantisedDiatomic.position_above_surface!(r, height.*direction)
    @test r == [0 0; height height; 0 0]
end

@testset "velocity_from_energy" begin
    v = QuantisedDiatomic.velocity_from_energy([100, 100], 10u"eV")
    v2 = QuantisedDiatomic.velocity_from_energy([100, 100], austrip(10u"eV"))
    @test v ≈ v2
end

@testset "apply_translational_impulse!" begin
    masses = [100, 100]
    v = rand(3, 2) ./ 100
    v_before = copy(v)
    e = 10u"eV"
    direction = [0, 0.5, 0.5]
    QuantisedDiatomic.apply_translational_impulse!(v, masses, e, direction)
    @test v_before[1,:] == v[1,:]
    @test all(v_before[2:3,:] .< v[2:3,:])

    νi, Ji = quantise_diatomic(sim, v, r)
    QuantisedDiatomic.apply_translational_impulse!(v, sim.atoms.masses, 1, rand(3))
    νf, Jf = quantise_diatomic(sim, v, r)
    @test νi ≈ νf
    @test Ji ≈ Jf
end
