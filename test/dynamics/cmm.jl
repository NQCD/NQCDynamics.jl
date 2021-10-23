
using Test
using NonadiabaticMolecularDynamics
using FiniteDiff

function test_motion!(sim::Simulation{<:CMM2}, u)
    f(x) = DynamicsUtils.classical_hamiltonian(sim, x)

    grad = FiniteDiff.finite_difference_gradient(f, u)

    du = zero(u)
    DynamicsMethods.motion!(du, u, sim, 0.0)

    @test DynamicsUtils.get_positions(du) ≈ DynamicsUtils.get_velocities(grad) ./ sim.atoms.masses' rtol=1e-3
    @test DynamicsUtils.get_velocities(du) ≈ -DynamicsUtils.get_positions(grad) ./ sim.atoms.masses' rtol=1e-3
    @test DynamicsMethods.MappingVariableMethods.get_mapping_positions(du) ≈ DynamicsMethods.MappingVariableMethods.get_mapping_momenta(grad) rtol=1e-3
    @test DynamicsMethods.MappingVariableMethods.get_mapping_momenta(du) ≈ -DynamicsMethods.MappingVariableMethods.get_mapping_positions(grad) rtol=1e-3
end

sim = Simulation{CMM2}(Atoms(1), DoubleWell())

v = randn(1,1)
r = randn(1,1)
u = DynamicsVariables(sim, v, r, 1)

test_motion!(sim, u)

sol = run_trajectory(u, (0, 100.0), sim; output=(:hamiltonian, :position), reltol=1e-10, abstol=1e-10)
@test sol.hamiltonian[1] ≈ sol.hamiltonian[end] rtol=1e-3
