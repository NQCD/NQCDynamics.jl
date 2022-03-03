using Test
using NQCDynamics
using Unitful

atoms = Atoms([:C, :H])
classical_sim = Simulation(atoms, Harmonic())
ring_polymer_sim = RingPolymerSimulation(atoms, Harmonic(), 10; temperature=100u"K")
sims = [classical_sim, ring_polymer_sim]

name(sim) = typeof(sim).name.wrapper

@testset "classical_potential_energy : $(name(sim))" for sim in sims
    u = DynamicsVariables(sim, randn(size(sim)), randn(size(sim)))
    r = DynamicsUtils.get_positions(u)
    @show DynamicsUtils.classical_potential_energy(sim, r)
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
