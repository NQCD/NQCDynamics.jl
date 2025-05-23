using Test
using NQCDynamics
using LinearAlgebra
using StaticArrays
using Statistics: mean
using StatsBase: sem
using DataFrames, CSV
using Interpolations
import JSON

benchmark_dir = get(ENV, "BENCHMARK_OUTPUT_DIR", "tmp/nqcd_benchmark")
benchmark_results = Dict{String, Any}("title_for_plotting" => "RPCME Tests")

ħω = 0.003
Γ = 0.003
kT = 0.03
g = 0.005
Ed = 0.0
atoms = Atoms(1 / ħω)

"""
Model used in J. Chem. Phys. 142, 084110 (2015).

Below we verify the implementation by reproducing the data in Fig. 2.
"""
struct TestModel{T} <: NQCModels.DiabaticModels.DiabaticModel
    ħω::T
    Ed::T
    g::T
    Γ::T
end

NQCModels.nstates(::TestModel) = 2
NQCModels.ndofs(::TestModel) = 1

function NQCModels.potential(model::TestModel, r::Real)
    (; ħω, Ed, g) = model
    potential = ħω * r^2 / 2
    V11 = potential
    V22 = potential + Ed + sqrt(2) * g * r
    V12 = sqrt(Γ / 2π)
    return Hermitian(SMatrix{2,2}(V11, V12, V12, V22))
end

function NQCModels.derivative(model::TestModel, r::Real)
    (; ħω, g) = model
    potential = ħω * r
    V11 = potential
    V22 = potential + sqrt(2) * g
    V12 = 0
    return Hermitian(SMatrix{2,2}(V11, V12, V12, V22))
end

model = TestModel(ħω, Ed, g, Γ)
n_beads = 4

@testset "Phonon relaxation Tᵢ=$(T)T" for T in [0.2, 2.0, 5.0]
    sim = RingPolymerSimulation{CME}(atoms, model, n_beads; temperature=kT)
    v = zeros(1, 1, n_beads)
    kTinitial = kT * T

    r = PositionHarmonicRingPolymer{Float64}(ħω, 1 / kTinitial, 1 / ħω, (1, 1, n_beads))
    v = VelocityBoltzmann(kTinitial * n_beads, atoms.masses[1])
    distribution = DynamicalDistribution(v, r, (1, 1, n_beads)) * PureState(1, Diabatic())

    dyn_test = @timed begin
        run_dynamics(sim, (0.0, 200 / Γ), distribution; trajectories=250,
            output=(OutputKineticEnergy, OutputTotalEnergy, OutputPotentialEnergy, OutputDiscreteState, OutputSpringEnergy, OutputCentroidKineticEnergy),
            abstol=1e-12, reltol=1e-12, saveat=2 / Γ, dt=1 / ħω / 10
        )
    end
    output = dyn_test.value
    
    avg = mean(o[:OutputCentroidKineticEnergy] for o in output) ./ kT
    err = zero(avg)
    for i in eachindex(err)
        err[i] = sem(o[:OutputCentroidKineticEnergy][i] for o in output; mean=avg[i] * kT) / kT
    end

    data = CSV.read(joinpath(@__DIR__, "reference_data", "$(T)T.csv"), DataFrame; header=false)
    itp = linear_interpolation(data[!, 1], data[!, 2]; extrapolation_bc=Line())
    for i in eachindex(avg)
        t = output[1][:Time][i] * Γ
        true_value = itp(t)
        @test isapprox(true_value, avg[i]; atol=10err[i], rtol=0.2)
    end
    # Uncomment to see the comparison if the tests start failing
    # using Plots
    # p = plot()
    # plot!(data[!,1], data[!,2]; label="$T")
    # plot!(output[1][:Time] .* Γ, avg, yerr=err)
    # plot!(output[1][:Time] .* Γ, itp.(output[1][:Time] .* Γ))
    # display(p)
    benchmark_results["Phonon relaxation"] = Dict("Time" => dyn_test.time, "Allocs" => dyn_test.bytes)
end

# Make benchmark directory if it doesn't already exist.
if !isdir(benchmark_dir)
    mkpath(benchmark_dir)
    @info "Benchmark data ouput directory created at $(benchmark_dir)."
else
    @info "Benchmark data ouput directory exists at $(benchmark_dir)."
end

# Output benchmarking dict
output_file = open("$(benchmark_dir)/rpcme.json", "w")
JSON.print(output_file, benchmark_results)
close(output_file)
