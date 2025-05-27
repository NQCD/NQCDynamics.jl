using Test
using NQCDynamics
using OrdinaryDiffEq: Tsit5
using Random: seed!
import JSON

benchmark_dir = get(ENV, "BENCHMARK_OUTPUT_DIR", "tmp/nqcd_benchmark")
benchmark_results = Dict{String, Any}("title_for_plotting" => "BCBwithTsit5 Tests")

@testset "Ehrenfest" begin
    sim = RingPolymerSimulation{Ehrenfest}(Atoms(2000), TullyModelOne(), 10; temperature=1e-3)
    u = DynamicsVariables(sim, fill(10/2000, size(sim)), fill(-5, size(sim)) .+ randn(size(sim)), PureState(1))
    dt = 0.1
    
    dyn_test = @timed run_dynamics(sim, (0, 2000.0), u; output=OutputDynamicsVariables, algorithm=DynamicsMethods.IntegrationAlgorithms.BCBwithTsit5(Tsit5()), saveat=0:10:2000, dt=dt)
    sol = dyn_test.value
    
    sol1 = run_dynamics(sim, (0, 2000.0), u; output=OutputDynamicsVariables, algorithm=Tsit5(), saveat=0:10:2000, abstol=1e-6, reltol=1e-6)

    @test sol[:OutputDynamicsVariables] ≈ sol1[:OutputDynamicsVariables] rtol=1e-3
    benchmark_results["Ehrenfest"] = Dict("Time" => dyn_test.time, "Allocs" => dyn_test.bytes)
end

@testset "FSSH" begin
    sim = RingPolymerSimulation{FSSH}(Atoms(2000), TullyModelOne(), 10; temperature=1e-3)
    u = DynamicsVariables(sim, fill(20/2000, size(sim)), fill(-5, size(sim)) .+ randn(size(sim)), PureState(1))
    dt = 0.1

    seed!(1)
    dyn_test = @timed run_dynamics(sim, (0, 2000.0), u; output=OutputDynamicsVariables, algorithm=DynamicsMethods.IntegrationAlgorithms.BCBwithTsit5(Tsit5()), saveat=0:10:2000, dt=dt)
    sol = dyn_test.value
    seed!(1)
    sol1 = run_dynamics(sim, (0, 2000.0), u; output=OutputDynamicsVariables, algorithm=Tsit5(), saveat=0:10:2000, dt=dt, adaptive=false)

    @test sol[:OutputDynamicsVariables] ≈ sol1[:OutputDynamicsVariables] rtol=1e-3
    benchmark_results["FSSH"] = Dict("Time" => dyn_test.time, "Allocs" => dyn_test.bytes)
end

# Make benchmark directory if it doesn't already exist.
if !isdir(benchmark_dir)
    mkpath(benchmark_dir)
    @info "Benchmark data ouput directory created at $(benchmark_dir)."
else
    @info "Benchmark data ouput directory exists at $(benchmark_dir)."
end

# Output benchmarking dict
output_file = open("$(benchmark_dir)/BCBwithTsit5.json", "w")
JSON.print(output_file, benchmark_results)
close(output_file)
