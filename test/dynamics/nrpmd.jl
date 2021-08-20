using Test
using NonadiabaticMolecularDynamics
using FiniteDiff
using Random
using OrdinaryDiffEq
using DiffEqDevTools
Random.seed!(1)

@test Dynamics.NRPMD{Float64}(10) isa Dynamics.NRPMD
atoms = Atoms(1.0)
sim = RingPolymerSimulation{NRPMD}(atoms, NonadiabaticModels.DoubleWell(), 10; DoFs=1, temperature=1e-1)

v = zeros(sim.DoFs, length(sim.atoms), length(sim.beads))
r = randn(sim.DoFs, length(sim.atoms), length(sim.beads))
u = DynamicsVariables(sim, v, r, 2)

@testset "get_diabatic_population" begin
    population = Dynamics.get_diabatic_population(sim, u)
    @test population[1] ≈ 0 atol=1e-10
    @test population[2] ≈ 1
end

@testset "get_adiabatic_population" begin
    sim = RingPolymerSimulation{NRPMD}(atoms, DoubleWell(), 10; DoFs=1, temperature=1e-1)
    u = DynamicsVariables(sim, zeros(1,1,10), zeros(1,1,10), 2)
    population = Dynamics.get_diabatic_population(sim, u)
    @test population[1] ≈ 0 atol=1e-10
    @test population[2] ≈ 1
    population = Dynamics.get_adiabatic_population(sim, u)
    @test population[1] ≈ 0.5
    @test population[2] ≈ 0.5
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

sol = Dynamics.run_trajectory(u, (0, 10.0), sim; output=(:hamiltonian, :position), dt=1e-2)
@test sol.hamiltonian[1] ≈ sol.hamiltonian[end] rtol=1e-2

@testset "MInt algorithm" begin
    tspan=(0, 20.0)
    prob = Dynamics.create_problem(u, tspan, sim)
    dts = (1/2) .^ (14:-1:8)
    setup = Dict(:alg => Feagin12(), :adaptive=>true, :reltol=>1e-14, :abstol=>1e-14)
    res = analyticless_test_convergence(dts, prob, Dynamics.MInt(), setup)
    @test res.𝒪est[:final] ≈ 2 atol=0.1
end
