using Test
using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.InitialConditions.QuantisedDiatomic
using LinearAlgebra: norm
using Unitful
using UnitfulAtomic

atoms = Atoms([:N, :O])
model = DiatomicHarmonic()
sim = Simulation(atoms, model)
surface = QuantisedDiatomic.SurfaceParameters(sim.atoms.masses, [1, 2], Matrix{Float64}(undef, 3, 0), 10.0, [0, 0, 1.0])

@testset "separate/combine_slab_and_molecule" begin
    atom_indices = [2, 3]
    r0 = [
        1 2 3 4 5;
        6 7 8 9 10;
        11 12 13 14 15
    ]
    molecule, slab = QuantisedDiatomic.separate_slab_and_molecule(atom_indices, r0)
    @test molecule == [2 3; 7 8; 12 13]
    @test slab == [1 4 5; 6 9 10; 11 14 15]

    r = QuantisedDiatomic.combine_slab_and_molecule(atom_indices, molecule, slab)
    @test r == r0

    atom_indices = [1, 2]
    r0 = [1 2; 3 4; 5 6]
    molecule, slab = QuantisedDiatomic.separate_slab_and_molecule(atom_indices, r0)
    @test molecule == r0

    r = QuantisedDiatomic.combine_slab_and_molecule(atom_indices, molecule, slab)
    @test r == r0
end

@testset "calculate_diatomic_energy" begin
    @test QuantisedDiatomic.calculate_diatomic_energy(model, 1.0, surface) ≈ 0.0 atol=1e-3
    @test QuantisedDiatomic.calculate_diatomic_energy(model, 2.0, surface) ≈ 0.5
end

@testset "calculate_force_constant" begin
    @test QuantisedDiatomic.calculate_force_constant(model, surface)[1] ≈ 1 atol=0.1
    model2 = NonadiabaticModels.DiatomicHarmonic(r₀=6.3)
    @test QuantisedDiatomic.calculate_force_constant(model2, surface)[1] ≈ 1 atol=0.1
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
    QuantisedDiatomic.position_above_surface!(r, height, InfiniteCell())
    @test r == [0 0; 0 0; height height]
end

@testset "velocity_from_energy" begin
    v = QuantisedDiatomic.velocity_from_energy([100, 100], 10u"eV")
    v2 = QuantisedDiatomic.velocity_from_energy([100, 100], austrip(10u"eV"))
    @test v ≈ v2
end

@testset "apply_translational_impulse!" begin
    masses = [100, 100]
    r = rand(3, 2)
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
