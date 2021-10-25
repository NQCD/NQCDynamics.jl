
using Test
using NonadiabaticMolecularDynamics
using FiniteDiff

function test_motion!(sim::Simulation{<:eCMM}, u)
    f(x) = DynamicsUtils.classical_hamiltonian(sim, x)

    grad = FiniteDiff.finite_difference_gradient(f, u)

    du = zero(u)
    DynamicsMethods.motion!(du, u, sim, 0.0)

    @test DynamicsUtils.get_positions(du) ≈ DynamicsUtils.get_velocities(grad) ./ sim.atoms.masses' rtol=1e-3
    @test DynamicsUtils.get_velocities(du) ≈ -DynamicsUtils.get_positions(grad) ./ sim.atoms.masses' rtol=1e-3
    @test DynamicsMethods.MappingVariableMethods.get_mapping_positions(du) ≈ DynamicsMethods.MappingVariableMethods.get_mapping_momenta(grad) rtol=1e-3
    @test DynamicsMethods.MappingVariableMethods.get_mapping_momenta(du) ≈ -DynamicsMethods.MappingVariableMethods.get_mapping_positions(grad) rtol=1e-3
end

sim = Simulation{eCMM}(Atoms(1), DoubleWell(); γ=0)
sim1 = Simulation{eCMM}(Atoms(1), DoubleWell(); γ=0.5)

v = randn(1,1)
r = randn(1,1)
u = DynamicsVariables(sim, v, r, SingleState(1))
u1 = DynamicsVariables(sim1, v, r, SingleState(1))

test_motion!(sim, u)
test_motion!(sim1, u1)

sol = run_trajectory(u, (0, 100.0), sim; output=(:hamiltonian, :position, :population), reltol=1e-10, abstol=1e-10)
@test sol.hamiltonian[1] ≈ sol.hamiltonian[end] rtol=1e-3
@test all(isapprox.(sum.(sol.population), 1, rtol=1e-3))

sol = run_trajectory(u1, (0, 100.0), sim1; output=(:hamiltonian, :position, :population), reltol=1e-10, abstol=1e-10)
@test sol.hamiltonian[1] ≈ sol.hamiltonian[end] rtol=1e-3
@test all(isapprox.(sum.(sol.population), 1, rtol=1e-3))
