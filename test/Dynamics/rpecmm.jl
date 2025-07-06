
using Test
using NQCDynamics
using FiniteDiff
using LinearAlgebra: norm
using OrdinaryDiffEq: Tsit5
import JSON

benchmark_dir = get(ENV, "BENCHMARK_OUTPUT_DIR", "tmp/nqcd_benchmark")
benchmark_results = Dict{String, Any}("title_for_plotting" => "RPeCMM Tests")

function test_motion!(sim::RingPolymerSimulation{<:eCMM}, u)
    f(x) = DynamicsUtils.classical_hamiltonian(sim, x)

    grad = FiniteDiff.finite_difference_gradient(f, u)

    du = zero(u)
    DynamicsMethods.motion!(du, u, sim, 0.0)

    @test DynamicsUtils.get_positions(du) ≈ DynamicsUtils.get_velocities(grad) ./ sim.atoms.masses' rtol=1e-3
    @test DynamicsUtils.get_velocities(du) ≈ -DynamicsUtils.get_positions(grad) ./ sim.atoms.masses' rtol=1e-3
    @test DynamicsMethods.MappingVariableMethods.get_mapping_positions(du) ≈ DynamicsMethods.MappingVariableMethods.get_mapping_momenta(grad) rtol=1e-3
    @test DynamicsMethods.MappingVariableMethods.get_mapping_momenta(du) ≈ -DynamicsMethods.MappingVariableMethods.get_mapping_positions(grad) rtol=1e-3
end

nbeads = 10
sim = RingPolymerSimulation{eCMM}(Atoms(1), DoubleWell(), nbeads; γ=0.0, temperature=0.01)
sim1 = RingPolymerSimulation{eCMM}(Atoms(1), DoubleWell(), nbeads; γ=0.5, temperature=0.01)

v = randn(1,1,nbeads)
r = zeros(1,1,nbeads)
u = DynamicsVariables(sim, v, r, PureState(1))
u1 = DynamicsVariables(sim1, v, r, PureState(1))

test_motion!(sim, u)
test_motion!(sim1, u1)

dyn_test = @timed run_dynamics(sim, (0, 100.0), u; output=(OutputTotalEnergy, OutputMappingPosition, OutputMappingMomentum), dt=0.01)
sol = dyn_test.value
@test sol[:OutputTotalEnergy][1] ≈ sol[:OutputTotalEnergy][end] rtol=1e-3
benchmark_results["RPeCMM γ=0.0"] = Dict("Time" => dyn_test.time, "Allocs" => dyn_test.bytes)
qmap = sol[:OutputMappingPosition]
pmap = sol[:OutputMappingMomentum]
total_population = sum.(DynamicsMethods.MappingVariableMethods.mapping_kernel.(qmap, pmap, sim.method.γ))
@test all(isapprox.(total_population, 1, rtol=1e-2))

dyn_test = @timed run_dynamics(sim1, (0, 100.0), u1; output=(OutputTotalEnergy, OutputMappingPosition, OutputMappingMomentum), dt=0.01)
sol = dyn_test.value
@test sol[:OutputTotalEnergy][1] ≈ sol[:OutputTotalEnergy][end] rtol=1e-3
benchmark_results["RPeCMM γ=0.5"] = Dict("Time" => dyn_test.time, "Allocs" => dyn_test.bytes)
qmap = sol[:OutputMappingPosition]
pmap = sol[:OutputMappingMomentum]
total_population = sum.(DynamicsMethods.MappingVariableMethods.mapping_kernel.(qmap, pmap, sim1.method.γ))
@test all(isapprox.(total_population, 1, rtol=1e-2))

@testset "Population correlation" begin
    # Tests that the initial state correlated with itself is 1 and correlated with other state is 0.
    gams = -0.5:0.5:1.5
    @testset "Gamma = $γ" for γ in gams[2:end]
        sim = RingPolymerSimulation{eCMM}(Atoms(1), DoubleWell(), 100; γ=γ)
        out = zeros(2, 2)
        n = 1e5
        correlation = TimeCorrelationFunctions.PopulationCorrelationFunction(sim, Diabatic())
        normalisation = TimeCorrelationFunctions.evaluate_normalisation(sim, correlation)

        for i=1:n
            u = DynamicsVariables(sim, 0, 0, PureState(1))
            K = Estimators.initial_diabatic_population(sim, u)
            Kinv = Estimators.diabatic_population(sim, u)
            out .+= normalisation * K * Kinv'
        end
        @test out ./ n ≈ [1 0; 0 1] atol=0.2
    end
end

@testset "Algorithm comparison" begin
    sol = run_dynamics(sim, (0, 10.0), u; output=OutputDynamicsVariables, dt=1e-2, algorithm=DynamicsMethods.IntegrationAlgorithms.RingPolymerMInt())
    sol1 = run_dynamics(sim, (0, 10.0), u; output=OutputDynamicsVariables, algorithm=Tsit5(), reltol=1e-10, abstol=1e-10, saveat=sol[:Time])
    @test sol[:OutputDynamicsVariables] ≈ sol1[:OutputDynamicsVariables] rtol=1e-2
    sol = run_dynamics(sim1, (0, 10.0), u; output=OutputDynamicsVariables, dt=1e-2, algorithm=DynamicsMethods.IntegrationAlgorithms.RingPolymerMInt())
    sol1 = run_dynamics(sim1, (0, 10.0), u; output=OutputDynamicsVariables, algorithm=Tsit5(), reltol=1e-10, abstol=1e-10, saveat=sol[:Time])
    @test sol[:OutputDynamicsVariables] ≈ sol1[:OutputDynamicsVariables] rtol=1e-2
end

# Make benchmark directory if it doesn't already exist.
if !isdir(benchmark_dir)
    mkpath(benchmark_dir)
    @info "Benchmark data ouput directory created at $(benchmark_dir)."
else
    @info "Benchmark data ouput directory exists at $(benchmark_dir)."
end

# Output benchmarking dict
output_file = open("$(benchmark_dir)/RPeCMM.json", "w")
JSON.print(output_file, benchmark_results)
close(output_file)
