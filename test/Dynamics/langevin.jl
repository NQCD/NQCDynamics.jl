using Test
using NQCDynamics
using JSON

benchmark_dir = get(ENV, "BENCHMARK_OUTPUT_DIR", "/tmp/nqcd_benchmark")
benchmark_results = Dict{String, Any}("title_for_plotting" => "Langevin")

atoms = Atoms([:H])
model = Harmonic()

sim = Simulation{Langevin}(atoms, model; temperature=1, γ=1)

z = DynamicsVariables(sim, randn(1,1), randn(1,1))

dyn_test = @timed run_dynamics(sim, (0.0, 500.0), z; output=OutputDynamicsVariables, dt=1)
solution = dyn_test.value
benchmark_results["Langevin"] = Dict("Time" => dyn_test.time, "Allocs" => dyn_test.bytes)

# Output benchmarking dict
output_file = open("$(benchmark_dir)/langevin.json", "w")
JSON.print(output_file, benchmark_results)
close(output_file)
