using Test
using NQCDynamics
using NQCCalculators
using Unitful

atoms = Atoms([:C, :H])
classical_sim = Simulation(atoms, Harmonic())
ring_polymer_sim = RingPolymerSimulation(atoms, Harmonic(), 10; temperature=100u"K")
sims = [classical_sim, ring_polymer_sim]

name(sim) = typeof(sim).name.wrapper

@testset "classical_potential_energy : $(name(sim))" for sim in sims
    u = DynamicsVariables(sim, randn(size(sim)), randn(size(sim)))
    r = DynamicsUtils.get_positions(u)
    @test DynamicsUtils.classical_potential_energy(sim, r) == DynamicsUtils.classical_potential_energy(sim, u)
end

@testset "classical_kinetic_energy : $(name(sim))" for sim in sims
    u = DynamicsVariables(sim, randn(size(sim)), randn(size(sim)))
    v = DynamicsUtils.get_velocities(u)
    @test DynamicsUtils.classical_kinetic_energy(sim, v) == DynamicsUtils.classical_kinetic_energy(sim, u)
end

@testset "divide_by_mass! : $(name(sim))" for sim in sims
    dv = randn(size(sim))
    true_value = copy(dv)
    for I in CartesianIndices(dv)
        true_value[I] /= masses(sim)[I[2]]
    end
    DynamicsUtils.divide_by_mass!(dv, masses(sim))
    @test true_value ≈ dv
end

rs = (fill(0.0, 1, 1), fill(0.0, 1, 1, 1))
caches = (
    NQCCalculators.Create_Cache(DoubleWell(), 1, Float64),
    NQCCalculators.Create_Cache(DoubleWell(), 1, 1, Float64),
    )
@testset "transform_density!" for (r, cache) ∈ zip(rs, caches)
    density = [1.0 0; 0 0]
    DynamicsUtils.transform_density!(density, cache, r, :to_adiabatic)
    @test density ≈ [0.5 0.5; 0.5 0.5]
    DynamicsUtils.transform_density!(density, cache, r, :to_diabatic)
    @test density ≈ [1.0 0; 0 0]
    @test_throws ArgumentError DynamicsUtils.transform_density!(density, cache, r, :blah)
end

@testset "initialise_adiabatic_density_matrix" for (r, cache) ∈ zip(rs, caches)
    @testset "Adiabatic" begin
        electronics = PureState(1, Adiabatic())
        density = DynamicsUtils.initialise_adiabatic_density_matrix(electronics, cache, r)
        @test density[1,1] == 1
    end

    @testset "Diabatic" begin
        electronics = PureState(1, Diabatic())
        density = DynamicsUtils.initialise_adiabatic_density_matrix(electronics, cache, r)
        @test all(density .≈ 0.5)
    end
end
