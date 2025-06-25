using NQCDynamics
using Test
using Unitful
using ComponentArrays
import JSON

include("utils.jl")

atoms = Atoms([:H])
model = NQCModels.Harmonic()

benchmark_dir = get(ENV, "BENCHMARK_OUTPUT_DIR", "tmp/nqcd_benchmark") # have made this director explicitly, if we want to be able to change this directory will need to put in a mkdir check
benchmark_results = Dict{String, Any}("title_for_plotting" => "Classical Tests")

@testset "Classical" begin
    sim = Simulation{Classical}(atoms, model)

    test_velocity!(sim)
    test_acceleration!(sim)
    test_motion!(sim)

    v = get_blank(sim)
    r = get_blank(sim)
    u0 = ComponentVector(v=v, r=r)

    dyn_test = @timed run_dynamics(sim, (0.0, 1000.0), u0; dt=0.1, output=(OutputTotalEnergy))
    sol = dyn_test.value
    @test isapprox(sol[:OutputTotalEnergy][1], sol[:OutputTotalEnergy][end], rtol=1e-2)
    benchmark_results["Classical"] = Dict("Time" => dyn_test.time, "Allocs" => dyn_test.bytes)
end

@testset "Ring polymer classical" begin
    sim = RingPolymerSimulation{Classical}(atoms, model, 10; temperature=10000u"K")

    test_velocity!(sim)
    test_acceleration!(sim)
    test_motion!(sim)

    v = get_blank(sim)
    r = get_blank(sim)
    u0 = ComponentVector(v=v, r=r)

    dyn_test = @timed run_dynamics(sim, (0.0, 1000.0), u0; dt=0.1, output=(OutputTotalEnergy))
    sol = dyn_test.value
    @test isapprox(sol[:OutputTotalEnergy][1], sol[:OutputTotalEnergy][end], rtol=1e-2)
    benchmark_results["Ring Polymer Classical"] = Dict("Time" => dyn_test.time, "Allocs" => dyn_test.bytes)
end

@testset "Fermion model classical adiabatic dynamics" begin
    model = AdiabaticStateSelector(WideBandBath(ErpenbeckThoss(;Γ=0.1u"eV"), bandmin=-5u"eV", bandmax=5u"eV", step=0.2u"eV"), 1)
    sim = Simulation{Classical}(atoms, model; temperature=300u"K")
    v = rand(VelocityBoltzmann(300u"K", atoms.masses, size(sim)))
    r = model.quantum_model.model.morse.x₀ |> hcat
    u0 = DynamicsVariables(sim, v, r)
    dyn_test = @timed run_dynamics(sim, (0.0, 900.0u"fs"), u0; dt=1u"fs", output=(OutputTotalEnergy))
    sol = dyn_test.value
    @test isapprox(sol[:OutputTotalEnergy][1], sol[:OutputTotalEnergy][end], atol=1e-2)
    benchmark_results["Fermion model classical adiabatic dynamics"] = Dict("Time" => dyn_test.time, "Allocs" => dyn_test.bytes)
end

@testset "Fermion model ring polymer adiabatic dynamics" begin
    full_model = WideBandBath(ErpenbeckThoss(;Γ=0.1u"eV"), bandmin=-5u"eV", bandmax=5u"eV", step=0.2u"eV")
    model = AdiabaticStateSelector(full_model, 1)
    n_beads = 10
    sim = RingPolymerSimulation{Classical}(atoms, model, n_beads; temperature=300u"K")
    v = VelocityBoltzmann(n_beads*300u"K", atoms.masses, size(sim)[1:2])
    r = model.quantum_model.model.morse.x₀
    d = DynamicalDistribution(v, r, size(sim))
    u0 = rand(d)
    dyn_test = @timed run_dynamics(sim, (0.0, 900.0u"fs"), u0; dt=1u"fs", output=(OutputTotalEnergy))
    sol = dyn_test.value
    @test sol[:OutputTotalEnergy][1] ≈ sol[:OutputTotalEnergy][end] rtol=1e-2
    benchmark_results["Fermion model ring polymer adiabatic dynamics"] = Dict("Time" => dyn_test.time, "Allocs" => dyn_test.bytes)
end

# Make benchmark directory if it doesn't already exist.
if !isdir(benchmark_dir)
    mkpath(benchmark_dir)
    @info "Benchmark data ouput directory created at $(benchmark_dir)."
else
    @info "Benchmark data ouput directory exists at $(benchmark_dir)."
end

# Output benchmarking dict
output_file = open("$(benchmark_dir)/classical.json", "w")
JSON.print(output_file, benchmark_results)
close(output_file)
