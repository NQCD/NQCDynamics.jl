using Test
using NQCDynamics
using Unitful
using UnitfulAtomic
using LinearAlgebra: diag
using ComponentArrays
using JSON

benchmark_dir = get(ENV, "BENCHMARK_OUTPUT_DIR", "tmp/nqcd_benchmark")
benchmark_results = Dict{String, Any}("title_for_plotting" => "MDEF Tests")

atoms = Atoms([:H, :H])
model = CompositeFrictionModel(Free(), ConstantFriction(1, atoms.masses[1]))
sim = Simulation{MDEF}(atoms, model; temperature=10u"K")

v = zeros(size(sim))
r = randn(size(sim))
u = ComponentVector(v=v, r=r)
du = zero(u)

@testset "friction!" begin
    gtmp = zeros(length(r), length(r))
    NQCDynamics.DynamicsMethods.ClassicalMethods.friction!(gtmp, r, sim, 0.0)
    @test all(diag(gtmp) .≈ 1.0)
end

dyn_test = @timed run_dynamics(sim, (0.0, 100.0), u; output=OutputDynamicsVariables, dt=1)
sol = dyn_test.value
@test sol[:OutputDynamicsVariables][1] ≈ u
benchmark_results["MDEF with ConstantFriction"] = Dict("Time" => dyn_test.time, "Allocs" => dyn_test.bytes)

f(t) = 100u"K"*exp(-ustrip(t))
model = CompositeFrictionModel(Harmonic(), RandomFriction(1))
sim = Simulation{MDEF}(atoms, model; temperature=f)

dyn_test = @timed run_dynamics(sim, (0.0, 100.0), u; output=OutputDynamicsVariables, dt=1)
sol = dyn_test.value
benchmark_results["MDEF with RandomFriction"] = Dict("Time" => dyn_test.time, "Allocs" => dyn_test.bytes)

# Make benchmark directory if it doesn't already exist.
if !isdir(benchmark_dir)
    mkpath(benchmark_dir)
    @info "Benchmark data ouput directory created at $(benchmark_dir)."
else
    @info "Benchmark data ouput directory exists at $(benchmark_dir)."
end

# Output benchmarking dict
output_file = open("$(benchmark_dir)/mdef.json", "w")
JSON.print(output_file, benchmark_results)
close(output_file)
