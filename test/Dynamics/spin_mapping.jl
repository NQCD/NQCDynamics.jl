using Test
using NQCDynamics
using FiniteDiff
using Random
using NQCDynamics: DynamicsMethods, DynamicsUtils
using NQCDynamics.DynamicsMethods.MappingVariableMethods
using OrdinaryDiffEq
using DiffEqDevTools
Random.seed!(1)

atoms = Atoms(1)
sim = Simulation{SpinMappingW}(atoms, DoubleWell())

v = randn(size(sim))
r = randn(size(sim))
u = DynamicsVariables(sim, v, r, PureState(2))

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
    sol = run_trajectory(u, (0, 10.0), sim; output=:hamiltonian, dt=1e-2)
    @test sol.hamiltonian[1] ≈ sol.hamiltonian[end] rtol=1e-2
end

@testset "Algorithm comparison" begin
    sol = run_trajectory(u, (0, 3.0), sim; dt=1e-2, algorithm=DynamicsMethods.IntegrationAlgorithms.MInt())
    sol1 = run_trajectory(u, (0, 3.0), sim; algorithm=Tsit5(), reltol=1e-10, abstol=1e-10, saveat=sol.t)
    @test sol.u ≈ sol1.u rtol=1e-2
end

@testset "MInt algorithm convergence" begin
    tspan=(0, 10.0)
    prob = DynamicsMethods.create_problem(u, tspan, sim)
    dts = 1 .// 2 .^(8:-1:2)

    alg = DynamicsMethods.IntegrationAlgorithms.MInt()
    test_alg = Vern9()
    setup = Dict(:alg => test_alg, :adaptive=>true, :abstol=>1e-8, :reltol=>1e-8)
    res = analyticless_test_convergence(dts, prob, alg, setup)
    @test res.𝒪est[:final] ≈ 2 atol=0.1
end
