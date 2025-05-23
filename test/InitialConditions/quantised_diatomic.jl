using Test
using NQCDynamics
using NQCDynamics.InitialConditions.QuantisedDiatomic
using LinearAlgebra: norm
using Unitful
using UnitfulAtomic

@testset "separate/combine_slab_and_molecule" begin

    @testset "3D, with slab" begin
        atom_indices = [2, 3]
        r0 = [
            1 2 3 4 5
            6 7 8 9 10
            11 12 13 14 15
        ]
        molecule, slab = QuantisedDiatomic.separate_slab_and_molecule(atom_indices, r0)
        @test molecule == [2 3; 7 8; 12 13]
        @test slab == [1 4 5; 6 9 10; 11 14 15]

        r = QuantisedDiatomic.combine_slab_and_molecule(atom_indices, molecule, slab)
        @test r == r0
    end

    @testset "3D, no slab" begin
        atom_indices = [1, 2]
        r0 = [1 2; 3 4; 5 6]
        molecule, slab = QuantisedDiatomic.separate_slab_and_molecule(atom_indices, r0)
        @test molecule == r0

        r = QuantisedDiatomic.combine_slab_and_molecule(atom_indices, molecule, slab)
        @test r == r0
    end

    @testset "1D, no slab" begin
        atom_indices = [1]
        r0 = [1;;]
        molecule, _ = QuantisedDiatomic.separate_slab_and_molecule(atom_indices, r0)
        @test molecule == r0
    end
end

@testset "calculate_force_constant" begin
    bond_lengths = 0.5:0.01:5.0
    model = NQCModels.DiatomicHarmonic(r₀=2.0)
    environment = QuantisedDiatomic.EvaluationEnvironment([1, 2], (3, 2), zeros(3, 0), 0.0, [0, 0, 1.0])
    binding_curve = QuantisedDiatomic.calculate_binding_curve(bond_lengths, model, environment)

    @test QuantisedDiatomic.calculate_force_constant(binding_curve) ≈ 1 atol = 0.1

    model = NQCModels.DiatomicHarmonic(r₀=3.0)
    binding_curve = QuantisedDiatomic.calculate_binding_curve(bond_lengths, model, environment)
    @test QuantisedDiatomic.calculate_force_constant(binding_curve) ≈ 1 atol = 0.1
end

@testset "subtract_centre_of_mass!" begin
    atoms = Atoms([:O, :O])
    r = rand(3, 2)
    r = QuantisedDiatomic.subtract_centre_of_mass(r, atoms.masses)
    @test r[:, 1] ≈ -r[:, 2]
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
    atoms = Atoms([:N, :O])
    model = DiatomicHarmonic(r₀=2)
    sim = Simulation(atoms, model)

    numbers = [2, 5]
    configs = generate_configurations(sim, numbers[1], numbers[2]; samples=10, translational_energy=1u"eV")
    for config in configs
        ν, J = quantise_diatomic(sim, config...)
        @test ν ≈ numbers[1] rtol = 1e-1
        @test J ≈ numbers[2] rtol = 1e-1
    end
end

@testset "position_above_surface!" begin
    r = zeros(3, 2)
    height = 10
    QuantisedDiatomic.position_above_surface!(r, [0, 0, height], InfiniteCell())
    @test r == [0 0; 0 0; height height]
end

@testset "velocity_from_energy" begin
    v = QuantisedDiatomic.velocity_from_energy([100, 100], 10u"eV")
    v2 = QuantisedDiatomic.velocity_from_energy([100, 100], austrip(10u"eV"))
    @test v ≈ v2
end

@testset "apply_translational_impulse!" begin
    atoms = Atoms([:N, :O])
    model = DiatomicHarmonic(r₀=2)
    sim = Simulation(atoms, model)

    r = [0.0 2.3; 0.0 0.0; 0.0 0.0]
    v = [1.0 2.0; 3.0 4.0; 5.0 6.0]
    v[:, 1] ./= sim.atoms.masses[1]
    v[:, 2] ./= sim.atoms.masses[2]
    v_before = copy(v)
    e = 10u"eV"
    direction = [0, 0.5, 0.5]
    QuantisedDiatomic.apply_translational_impulse!(v, sim.atoms.masses, e, direction)
    @test v_before[1, :] == v[1, :]
    @test all(v_before[2:3, :] .< v[2:3, :])

    νi, Ji = quantise_diatomic(sim, v, r; bond_lengths=0.1:0.01:10.0)
    QuantisedDiatomic.apply_translational_impulse!(v, sim.atoms.masses, 1, rand(3))
    νf, Jf = quantise_diatomic(sim, v, r; bond_lengths=0.1:0.01:10.0)
    @test νi ≈ νf
    @test Ji ≈ Jf
end

@testset "find_total_energy" begin
    environment = QuantisedDiatomic.EvaluationEnvironment([1, 2], (3, 2), zeros(3, 0), 0.0, [0, 0, 1.0])
    model = DiatomicHarmonic(r₀=2.0)
    bond_lengths = 0.5:0.01:5.0
    binding_curve = QuantisedDiatomic.calculate_binding_curve(bond_lengths, model, environment)
    μ = 1.0
    J = 2
    ν = 1
    V = QuantisedDiatomic.EffectivePotential(μ, J, binding_curve, true)
    E, bounds = QuantisedDiatomic.find_total_energy(V, ν)
end

@testset "find_integral_bounds" begin
    environment = QuantisedDiatomic.EvaluationEnvironment([1, 2], (3, 2), zeros(3, 0), 0.0, [0, 0, 1.0])
    model = DiatomicHarmonic(r₀=2)
    bond_lengths = 0.5:0.01:5.0
    binding_curve = QuantisedDiatomic.calculate_binding_curve(bond_lengths, model, environment)
    total_energy = 1
    J = 0
    μ = 1.0
    equilibrium_bond_length = 2
    V = QuantisedDiatomic.EffectivePotential(μ, J, binding_curve, true)
    r₁, r₂ = QuantisedDiatomic.find_integral_bounds(total_energy, V)
    @test r₁ < equilibrium_bond_length
    @test r₂ > equilibrium_bond_length
end

#= @testset "calculate_binding_curve" begin
    model = Morse()
    bond_lengths = 0.5:0.01:5.0
    environment = QuantisedDiatomic.EvaluationEnvironment([1], (1, 1), zeros(1, 0), 0.0, [0.0])
    binding_curve = QuantisedDiatomic.calculate_binding_curve(bond_lengths, model, environment)
    @test binding_curve.equilibrium_bond_length ≈ model.x₀
end =#

@testset "calculate_diatomic_energy" begin
    @testset "3D, no slab" begin
        model = DiatomicHarmonic(r₀=2)
        environment = QuantisedDiatomic.EvaluationEnvironment([1, 2], (3, 2), zeros(3, 0), 0.0, [0, 0, 1.0])
        @test QuantisedDiatomic.calculate_diatomic_energy(2.0, model, environment) ≈ 0.0 atol = 1e-3
        @test QuantisedDiatomic.calculate_diatomic_energy(3.0, model, environment) ≈ 0.5
    end

    #= @testset "1D, no slab" begin
        model = Harmonic(dofs=1)
        bond_length = 1.5
        environment = QuantisedDiatomic.EvaluationEnvironment([1], (1, 1), zeros(1, 0), 0.0, [0.0])
        @test QuantisedDiatomic.calculate_diatomic_energy(bond_length, model, environment) ≈ bond_length^2 / 2
    end =#
end

@testset "assemble_evaluation_geometry" begin
    bond_length = 4.3

    @testset "3D, no slab" begin
        environment = QuantisedDiatomic.EvaluationEnvironment([1, 2], (3, 4), rand(3, 2), 2.0, [0, 0, 1.0])
        r = QuantisedDiatomic.assemble_evaluation_geometry(bond_length, environment)
        @test norm(r[:, 1] .- r[:, 2]) == 4.3
    end

    #= @testset "1D, no slab" begin
        environment = QuantisedDiatomic.EvaluationEnvironment([1], (1, 1), zeros(1, 0), 0.0, [0.0])
        r = QuantisedDiatomic.assemble_evaluation_geometry(bond_length, environment)
        @test r == [4.3;;]
    end =#
end

@testset "build_molecule" begin
    bond_length = 4.3

    @testset "3D, no slab" begin
        environment = QuantisedDiatomic.EvaluationEnvironment([1, 2], (3, 4), rand(3, 2), 2.0, [0, 0, 1.0])
        r = QuantisedDiatomic.build_molecule(bond_length, environment)
        @test norm(r[:, 1] .- r[:, 2]) == 4.3
    end

    #= @testset "1D, no slab" begin
        environment = QuantisedDiatomic.EvaluationEnvironment([1], (1, 1), zeros(0, 0), 0.0, [0.0])
        r = QuantisedDiatomic.build_molecule(bond_length, environment)
        @test r == [4.3;;]
    end =#
end

@testset "fit_binding_curve and plot_binding_curve" begin
    model = Morse()
    bond_lengths = 0.5:0.01:5.0
    binding_curve = potential.(model, hcat.(bond_lengths))
    fit = QuantisedDiatomic.fit_binding_curve(bond_lengths, binding_curve)
    QuantisedDiatomic.plot_binding_curve(bond_lengths, binding_curve, fit)
end

@testset "generate_configurations 1D" begin
    sim = Simulation(Atoms(2000), Morse())

    @testset "forward and backward: ν = $quantum_number" for quantum_number in (0, 1, 2)
        r, v = QuantisedDiatomic.generate_1D_vibrations(Morse(), 100.0, quantum_number; samples=100)
        ν = QuantisedDiatomic.quantise_1D_vibration.(Morse(), 100.0, r, v)
        @test all(ν .== quantum_number)
    end

end
