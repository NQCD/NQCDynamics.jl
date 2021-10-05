using Test
using NonadiabaticMolecularDynamics
using FiniteDiff
using Random
using OrdinaryDiffEq
using DiffEqDevTools
using NonadiabaticMolecularDynamics: DynamicsMethods, DynamicsUtils
using NonadiabaticMolecularDynamics.DynamicsMethods.MappingVariableMethods
Random.seed!(1)

@test MappingVariableMethods.NRPMD{Float64}(10) isa MappingVariableMethods.NRPMD
atoms = Atoms(1.0)
sim = RingPolymerSimulation{NRPMD}(atoms, NonadiabaticModels.DoubleWell(), 10; temperature=1e-1)

v = zeros(size(sim))
r = randn(size(sim))
u = DynamicsVariables(sim, v, r, 2)

@testset "get_population" begin
    population = Estimators.diabatic_population(sim, u)
    @test population[1] ≈ 0 atol=1e-10
    @test population[2] ≈ 1
end

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

sol = run_trajectory(u, (0, 10.0), sim; output=(:hamiltonian, :position), dt=1e-2)
@test sol.hamiltonian[1] ≈ sol.hamiltonian[end] rtol=1e-2

@testset "MInt algorithm" begin
    tspan=(0, 20.0)
    prob = DynamicsMethods.create_problem(u, tspan, sim)
    dts = (1/2) .^ (14:-1:8)
    setup = Dict(:alg => Feagin12(), :adaptive=>true, :reltol=>1e-14, :abstol=>1e-14)
    res = analyticless_test_convergence(dts, prob, DynamicsMethods.IntegrationAlgorithms.MInt(), setup)
    @test res.𝒪est[:final] ≈ 2 atol=0.1
end
