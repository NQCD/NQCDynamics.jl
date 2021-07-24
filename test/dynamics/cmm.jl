
using Test
using NonadiabaticMolecularDynamics
using FiniteDiff

function test_motion!(sim::Simulation{<:CMM2}, u::Dynamics.MeyerMillerMappingVariables)
    f(x) = evaluate_hamiltonian(sim, x)

    grad = FiniteDiff.finite_difference_gradient(f, u)

    du = zero(u)
    Dynamics.motion!(du, u, sim, 0.0)

    @test get_positions(du) ≈ get_velocities(grad) ./ sim.atoms.masses' rtol=1e-3
    @test get_velocities(du) ≈ -get_positions(grad) ./ sim.atoms.masses' rtol=1e-3
    @test Dynamics.get_mapping_positions(du) ≈ Dynamics.get_mapping_momenta(grad) rtol=1e-3
    @test Dynamics.get_mapping_momenta(du) ≈ -Dynamics.get_mapping_positions(grad) rtol=1e-3
end

test_motion!(sim, u)

v = randn(1,1)
r = randn(1,1)
u = Dynamics.MeyerMillerMappingVariables(v, r, 2, 1)

sim = Simulation{CMM2}(Atoms(1), DoubleWell(); DoFs=1)

sol = Dynamics.run_trajectory(u, (0, 100.0), sim; output=(:hamiltonian, :position), reltol=1e-6)
@test sol.hamiltonian[1] ≈ sol.hamiltonian[end] rtol=1e-3
