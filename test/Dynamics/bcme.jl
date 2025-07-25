using Test
using NQCDynamics
using LinearAlgebra
using StaticArrays
using Distributions: Normal
using Statistics: mean
using StatsBase: sem
using DataFrames, CSV
using Interpolations
import JSON

benchmark_dir = get(ENV, "BENCHMARK_OUTPUT_DIR", "tmp/nqcd_benchmark")
benchmark_results = Dict{String, Any}("title_for_plotting" => "BCME Tests")

ħω = 0.003
Γ = 0.1
kT = 0.01
g = 0.015
Ed = g^2 / ħω
atoms = Atoms(1/ħω)
x1 = -sqrt(2) * g / ħω

"""
Model used in J. Chem. Phys. 142, 084110 (2015).

Below we verify the implementation by reproducing the data in Fig. 2.
"""
struct TestModel{T} <: NQCModels.QuantumModels.QuantumModel
    ħω::T
    Ed::T
    g::T
    Γ::T
end

NQCModels.nstates(::TestModel) = 2
NQCModels.ndofs(::TestModel) = 1

function NQCModels.potential!(model::TestModel, V::Hermitian, r::AbstractMatrix)
    (;ħω, Ed, g) = model
    potential = ħω*first(r)^2/2
    V11 = potential
    V22 = potential + Ed + sqrt(2)*g*first(r)
    V12 = sqrt(Γ/2π)
    V.data .= [V11 V12; V12 V22]
    return nothing
end

function NQCModels.derivative!(model::TestModel, D::AbstractMatrix{<:Hermitian}, r::AbstractMatrix)
    (;ħω, g) = model
    potential = ħω*first(r)
    V11 = potential
    V22 = potential + sqrt(2)*g
    V12 = 0
    D[1].data .= [V11 V12; V12 V22]
    return nothing
end

model = TestModel(ħω, Ed, g, Γ)

@testset "BCME" begin
    sim = Simulation{BCME}(atoms, model; temperature=kT, bandwidth=100.)
    βω = ħω / kT
    σ = sqrt(1 / βω)
    r = hcat(Normal(x1, σ))

    v = VelocityBoltzmann(kT, atoms.masses, (1,1))
    distribution = DynamicalDistribution(v, r, (1,1)) * PureState(1, Diabatic())

    dyn_test = @timed run_dynamics(sim, (0.0, 10000.0), distribution; trajectories=100, output=(OutputPosition), abstol=1e-12, reltol=1e-12, saveat=100)
    output = dyn_test.value
    
    avg = zeros(101)
    for i in eachindex(output)
        for j in eachindex(avg)
            avg[j] += output[i][:OutputPosition][j][1] / x1
        end
    end
    avg ./= length(output)
    err = zero(avg)
    for i in eachindex(err)
        err[i] = sem(o[:OutputPosition][i][1] / x1 for o in output; mean=avg[i])
    end

    data = CSV.read(joinpath(@__DIR__, "reference_data", "bcme.csv"), DataFrame; header=false)
    itp = linear_interpolation(data[!,1], data[!,2]; extrapolation_bc=Line())
    for i in eachindex(avg)
        t = output[1][:Time][i]
        true_value = itp(t)
        @test isapprox(true_value, avg[i]; atol=10err[i], rtol=0.1)
    end
    # Uncomment to see the comparison if the tests start failing
    # using Plots
    # p = plot()
    # plot!(data[!,1], data[!,2])
    # plot!(output[1][:Time], avg; yerr=err)
    # plot!(output[1][:Time], itp.(output[1][:Time]))
    # xlims!(0, 1e4)
    # display(p)
    benchmark_results["BCME"] = Dict("Time" => dyn_test.time, "Allocs" => dyn_test.bytes)
end

# Make benchmark directory if it doesn't already exist.
if !isdir(benchmark_dir)
    mkpath(benchmark_dir)
    @info "Benchmark data ouput directory created at $(benchmark_dir)."
else
    @info "Benchmark data ouput directory exists at $(benchmark_dir)."
end

# Output benchmarking dict
output_file = open("$(benchmark_dir)/BCME.json", "w")
JSON.print(output_file, benchmark_results)
close(output_file)
