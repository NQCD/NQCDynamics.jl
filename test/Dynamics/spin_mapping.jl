using Test
using NQCDynamics
using FiniteDiff
using Random
using NQCDynamics: DynamicsMethods, DynamicsUtils
using NQCDynamics.DynamicsMethods.MappingVariableMethods
# Random.seed!(1)

atoms = Atoms(1)
sim = Simulation{SpinMappingW}(atoms, DoubleWell())

v = randn(size(sim))
r = randn(size(sim))
u = DynamicsVariables(sim, v, r, SingleState(2))

@testset "motion! obeys Hamilton's equations" begin
    function test_motion!(sim, u)
        f(x) = DynamicsUtils.classical_hamiltonian(sim, x)

        grad = FiniteDiff.finite_difference_gradient(f, u)

        du = zero(u)
        DynamicsMethods.motion!(du, u, sim, 0.0)

        @test DynamicsUtils.get_positions(du) ≈ DynamicsUtils.get_velocities(grad) ./ sim.atoms.masses' rtol=1e-3
        @test DynamicsUtils.get_velocities(du) ≈ -DynamicsUtils.get_positions(grad) ./ sim.atoms.masses' rtol=1e-3
        @test MappingVariableMethods.get_mapping_positions(du) ≈ MappingVariableMethods.get_mapping_momenta(grad) rtol=1e-3
        @test MappingVariableMethods.get_mapping_momenta(du) ≈ -MappingVariableMethods.get_mapping_positions(grad) rtol=1e-3
    end

    test_motion!(sim, u)
end

@testset "Energy conservation" begin
    sol = run_trajectory(u, (0, 10.0), sim; output=:hamiltonian, atol=1e-10, reltol=1e-10)
    @test sol.hamiltonian[1] ≈ sol.hamiltonian[end] rtol=1e-2
end
