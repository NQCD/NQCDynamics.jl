using Test
using NQCDynamics
using LinearAlgebra
using StaticArrays
using Distributions: Normal
using Statistics: mean
using StatsBase: sem
using DataFrames, CSV
using Interpolations

ħω = 0.003
Γ = 0.003
kT = 0.03
g = 0.005
Ed = 0.0
atoms = Atoms(1/ħω)

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
    (;ħω, Ed, g) = model
    potential = ħω*r^2/2
    V11 = potential
    V22 = potential + Ed + sqrt(2)*g*r
    V12 = sqrt(Γ/2π)
    return Hermitian(SMatrix{2,2}(V11, V12, V12, V22))
end

function NQCModels.derivative(model::TestModel, r::Real)
    (;ħω, g) = model
    potential = ħω*r
    V11 = potential
    V22 = potential + sqrt(2)*g
    V12 = 0
    return Hermitian(SMatrix{2,2}(V11, V12, V12, V22))
end

model = TestModel(ħω, Ed, g, Γ)

@testset "Phonon relaxation Tᵢ=$(T)T" for T in [0.2, 2.0, 5.0]
    sim = Simulation{CME}(atoms, model; temperature=kT)
    v = zeros(1,1)
    r = zeros(1,1)
    kTinitial = kT * T
    βω = ħω / kTinitial
    σ = sqrt(1 / βω)
    r = Normal(0.0, σ)

    v = VelocityBoltzmann(kTinitial, atoms.masses[1])
    distribution = DynamicalDistribution(v, r, (1,1)) * PureState(1, Diabatic())

    output = run_dynamics(sim, (0.0, 200 / Γ), distribution; trajectories=500, output=(OutputKineticEnergy, OutputTotalEnergy, OutputDiscreteState), abstol=1e-8, reltol=1e-8, saveat=2/Γ)
    avg = mean(o[:OutputKineticEnergy] for o in output) ./ kT
    err = zero(avg)
    for i in eachindex(err)
        err[i] = sem(o[:OutputKineticEnergy][i] for o in output; mean=avg[i]*kT) / kT
    end

    data = CSV.read(joinpath(@__DIR__, "reference_data", "$(T)T.csv"), DataFrame; header=false)
    itp = linear_interpolation(data[!,1], data[!,2]; extrapolation_bc=Line())
    for i in eachindex(avg)
        t = output[1][:Time][i] * Γ
        true_value = itp(t)
        @test isapprox(true_value, avg[i]; atol=4err[i], rtol=0.2)
    end
    # Uncomment to see the comparison if the tests start failing
    # using Plots
    # p = plot()
    # plot!(data[!,1], data[!,2]; label="$T")
    # plot!(output[1][:Time] .* Γ, avg, yerr=err)
    # plot!(output[1][:Time] .* Γ, itp.(output[1][:Time] .* Γ))
    # display(p)
end

