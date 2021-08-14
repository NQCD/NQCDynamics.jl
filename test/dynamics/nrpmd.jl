using Test
using NonadiabaticMolecularDynamics
using FiniteDiff
using Random
Random.seed!(1)

@test Dynamics.NRPMD{Float64}(10) isa Dynamics.NRPMD
atoms = Atoms(1.0)
sim = RingPolymerSimulation{NRPMD}(atoms, NonadiabaticModels.DoubleWell(), 10; DoFs=1, temperature=1e-1)

v = zeros(sim.DoFs, length(sim.atoms), length(sim.beads))
r = randn(sim.DoFs, length(sim.atoms), length(sim.beads))
u = DynamicsVariables(sim, v, r, 2)

@testset "get_population" begin
    population = Dynamics.get_diabatic_population(sim, u)
    @test population[1] ≈ 0 atol=1e-10
    @test population[2] ≈ 1
end

function test_motion!(sim, u)
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

sol = Dynamics.run_trajectory(u, (0, 100.0), sim; output=(:hamiltonian, :position),
    reltol=1e-8, abstol=1e-8)
@test sol.hamiltonian[1] ≈ sol.hamiltonian[end] rtol=1e-2
