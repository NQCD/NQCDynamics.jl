
using Test
using NQCDynamics
using FiniteDiff
using LinearAlgebra: norm

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

sim = Simulation{eCMM}(Atoms(1), DoubleWell(); γ=0.0)
sim1 = Simulation{eCMM}(Atoms(1), DoubleWell(); γ=0.5)

v = randn(1,1)
r = randn(1,1)
u = DynamicsVariables(sim, v, r, PureState(1))
u1 = DynamicsVariables(sim1, v, r, PureState(1))

test_motion!(sim, u)
test_motion!(sim1, u1)

sol = run_dynamics(sim, (0, 100.0), u; output=(OutputTotalEnergy, OutputPosition, OutputDynamicsVariables), reltol=1e-10, abstol=1e-10)
@test sol[:OutputTotalEnergy][1] ≈ sol[:OutputTotalEnergy][end] rtol=1e-3
qmap = [u.qmap for u in sol[:OutputDynamicsVariables]]
pmap = [u.pmap for u in sol[:OutputDynamicsVariables]]
total_population = sum.(DynamicsMethods.MappingVariableMethods.mapping_kernel.(qmap, pmap, sim.method.γ))
@test all(isapprox.(total_population, 1, rtol=1e-3))

sol = run_dynamics(sim1, (0, 100.0), u1; output=(OutputTotalEnergy, OutputPosition, OutputDynamicsVariables), reltol=1e-10, abstol=1e-10)
@test sol[:OutputTotalEnergy][1] ≈ sol[:OutputTotalEnergy][end] rtol=1e-3
qmap = [u.qmap for u in sol[:OutputDynamicsVariables]]
pmap = [u.pmap for u in sol[:OutputDynamicsVariables]]
total_population = sum.(DynamicsMethods.MappingVariableMethods.mapping_kernel.(qmap, pmap, sim1.method.γ))
@test all(isapprox.(total_population, 1, rtol=1e-3))

@testset "generate_random_points_on_nsphere" begin
    points = DynamicsMethods.MappingVariableMethods.generate_random_points_on_nsphere(10, 1)
    @test norm(points) ≈ 1
    points = DynamicsMethods.MappingVariableMethods.generate_random_points_on_nsphere(10, 10)
    @test norm(points) ≈ 10
end

@testset "Population correlation" begin
    # Tests that the initial state correlated with itself is 1 and correlated with other state is 0.
    gams = -0.5:0.5:1.5
    for γ in gams[2:end]
        @testset "Gamma = $γ" begin
            sim = Simulation{eCMM}(Atoms(1), DoubleWell(); γ=γ)
            out = zeros(2, 2)
            n = 1e4
            correlation = TimeCorrelationFunctions.PopulationCorrelationFunction(sim, Diabatic())
            normalisation = TimeCorrelationFunctions.evaluate_normalisation(sim, correlation)

            for i=1:n
                u = DynamicsVariables(sim, 0, 0, PureState(1))
                K = Estimators.initial_diabatic_population(sim, u)
                Kinv = Estimators.diabatic_population(sim, u)
                out .+= normalisation * K * Kinv'
            end
            @test out ./ n ≈ [1 0; 0 1] atol=0.1
        end
    end
end

