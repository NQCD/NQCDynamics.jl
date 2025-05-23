using Test
using NQCDynamics
using Unitful
using UnitfulAtomic
using RecursiveArrayTools
using LinearAlgebra: diag
using StatsBase
using StochasticDiffEq
using ComponentArrays
using NQCDynamics: DynamicsMethods, DynamicsUtils
using NQCDynamics.DynamicsMethods: ClassicalMethods, IntegrationAlgorithms
import JSON

benchmark_dir = get(ENV, "BENCHMARK_OUTPUT_DIR", "tmp/nqcd_benchmark")
benchmark_results = Dict{String, Any}("title_for_plotting" => "RPMDEF Tests")

atoms = Atoms([:H, :C])
model = CompositeFrictionModel(Free(3), ConstantFriction(3, 1))
sim = RingPolymerSimulation{MDEF}(atoms, model, 3; temperature=10u"K")

v = RingPolymerArray(randn(size(sim)))
r = RingPolymerArray(randn(size(sim)))

@testset "friction!" begin
    gtmp = zeros(length(r), length(r))
    gtmp = zeros(ndofs(sim)*natoms(sim),ndofs(sim)*natoms(sim),nbeads(sim))
    F = ClassicalMethods.friction!(gtmp, r, sim, 0.0)
    hmass = sim.calculator.model.friction_model.γ/atoms.masses[1]
    cmass = sim.calculator.model.friction_model.γ/atoms.masses[2]
    for i=1:length(sim.beads)
        @test diag(F[:,:,i]) ≈ [hmass, hmass, hmass, cmass, cmass, cmass]
    end
end

prob = DynamicsMethods.create_problem(ComponentVector(v=v,r=r), (0.0, 0.5), sim)

@testset "step_C!" begin
    dt = 0.5
    c = RingPolymers.cayley_propagator(sim.beads, dt; half=true)
    ω_k = RingPolymers.get_matsubara_frequencies(length(sim.beads), sim.beads.ω_n)
    vbefore = copy(v)
    rbefore = copy(r)
    IntegrationAlgorithms.step_C!(v,r,c)
    for I in CartesianIndices(r)
        A = [0 1; -ω_k[I[3]] 0]
        a = exp(A*dt/2) * [rbefore[I], vbefore[I]]
        @test a ≈ [r[I], v[I]] rtol=1e-4
    end
end

@testset "step_O!" begin
    cache = StochasticDiffEq.alg_cache(IntegrationAlgorithms.BCOCB(), prob, ArrayPartition(v,r),0,0,sim,0,0,0,Float64,0,0,0,0,0,0,Val{true})
    integrator = init(prob,IntegrationAlgorithms.BCOCB();dt=0.5)
    IntegrationAlgorithms.step_O!(cache.friction, integrator, v, r, 0.0)
end

@testset "ThermalLangevin" begin
    atoms = Atoms([:H, :C])
    sim = RingPolymerSimulation{ThermalLangevin}(atoms, NQCModels.Free(), 3; temperature=100u"K", γ=0.01)

    v = RingPolymerArray(zeros(size(sim)))
    r = RingPolymerArray(zeros(size(sim)))

    dyn_test = @timed run_dynamics(sim, (0.0, 1e4), ArrayPartition(v,r); dt=1, output=OutputKineticEnergy)
    sol = dyn_test.value
    benchmark_results["ThermalLangevin"] = Dict("Time" => dyn_test.time, "Allocs" => dyn_test.bytes)
end

@testset "MDEF" begin
    atoms = Atoms([:H, :C])
    model = CompositeFrictionModel(Harmonic(), RandomFriction(1))
    sim = RingPolymerSimulation{MDEF}(atoms, model, 3; temperature=100u"K")

    v = RingPolymerArray(zeros(size(sim)))
    r = RingPolymerArray(zeros(size(sim)))

    dyn_test = @timed run_dynamics(sim, (0.0, 1e2), ArrayPartition(v,r); dt=1, output=OutputKineticEnergy)
    sol = dyn_test.value
    benchmark_results["RPMDEF"] = Dict("Time" => dyn_test.time, "Allocs" => dyn_test.bytes)
end

# Make benchmark directory if it doesn't already exist.
if !isdir(benchmark_dir)
    mkpath(benchmark_dir)
    @Info "Benchmark data ouput directory created at $(benchmark_dir)."
else
    @Info "Benchmark data ouput directory exists at $(benchmark_dir)."
end

# Output benchmarking dict
output_file = open("$(benchmark_dir)/RPMDEF.json", "w")
JSON.print(output_file, benchmark_results)
close(output_file)
