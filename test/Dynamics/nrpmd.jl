using Test
using NQCDynamics
using FiniteDiff
using Random
using OrdinaryDiffEq
using DiffEqDevTools
using NQCDynamics: DynamicsMethods, DynamicsUtils
using NQCDynamics.DynamicsMethods.MappingVariableMethods
Random.seed!(1)

@test MappingVariableMethods.NRPMD{Float64}(10, 0.0) isa MappingVariableMethods.NRPMD
atoms = Atoms(1)
sim = RingPolymerSimulation{NRPMD}(atoms, DoubleWell(), 10; temperature=1e-3)

v = randn(size(sim))
r = randn(size(sim))
u = DynamicsVariables(sim, v, r, PureState(2))

@testset "Population correlation" begin
    K = NQCModels.nstates(sim)
    out = zeros(K, K)
    n = 1e4
    correlation = TimeCorrelationFunctions.PopulationCorrelationFunction(sim, Diabatic())
    normalisation = TimeCorrelationFunctions.evaluate_normalisation(sim, correlation)

    for i=1:n
        u = DynamicsVariables(sim, v, r, PureState(1))
        K = Estimators.initial_diabatic_population(sim, u)
        Kinv = Estimators.diabatic_population(sim, u)
        out += normalisation * K * Kinv'
    end
   @test out ./ n ≈ [1 0; 0 0] atol=0.1
end

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

algs = (DynamicsMethods.IntegrationAlgorithms.RingPolymerMInt(), Tsit5())
@testset "Energy conservation $alg" for alg in algs
    sol = run_trajectory(u, (0, 10.0), sim; output=(:hamiltonian, :population), dt=1e-2, algorithm=alg, abstol=1e-8, reltol=1e-8)
    @test sol.hamiltonian[1] ≈ sol.hamiltonian[end] rtol=1e-2
end

@testset "Algorithm comparison" begin
    sol = run_trajectory(u, (0, 10.0), sim; dt=1e-2, algorithm=DynamicsMethods.IntegrationAlgorithms.RingPolymerMInt())
    sol1 = run_trajectory(u, (0, 10.0), sim; algorithm=Tsit5(), reltol=1e-10, abstol=1e-10, saveat=sol.t)
    @test sol.u ≈ sol1.u rtol=1e-2
end

@testset "MInt algorithm convergence" begin
    tspan=(0, 10.0)
    prob = DynamicsMethods.create_problem(u, tspan, sim)
    dts = 1 .// 2 .^(8:-1:4)

    alg = DynamicsMethods.IntegrationAlgorithms.RingPolymerMInt()
    test_alg = Vern9()
    setup = Dict(:alg => test_alg, :adaptive=>true, :abstol=>1e-8, :reltol=>1e-8)
    res = analyticless_test_convergence(dts, prob, alg, setup)
    @test res.𝒪est[:final] ≈ 2 atol=0.1
end
