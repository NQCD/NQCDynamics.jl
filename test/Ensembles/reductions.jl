
using NQCDynamics
using Test
using OrdinaryDiffEq
using HDF5

@testset "$reduction" for reduction in (MeanReduction(), SumReduction())
    prob = ODEProblem((u,p,t) -> 1.01u, [0.5], (0.0,1.0))
    prob_func(prob, i, repeat) = remake(prob,u0=[rand()])
    output = Ensembles.EnsembleSaver((OutputFinal,), true)
    ensemble = EnsembleProblem(prob, prob_func=prob_func, output_func=output, reduction=reduction)
    sol = solve(ensemble, Tsit5(), trajectories=1e3, batch_size=20, saveat=0.5)
    if reduction isa Ensembles.SumReduction
        @test sol.u[:OutputFinal][1] / 1e3 ≈ 0.5 * exp(1.01) rtol=1e-1
    else
        @test sol.u[:OutputFinal][1] ≈ 0.5 * exp(1.01) rtol=1e-1
    end
end

@testset "FileReduction" begin
    rm("test.h5", force=true)
    reduction = FileReduction("test.h5")

    output = (OutputPosition, OutputTotalEnergy)
    sim = Simulation(Atoms(1), Harmonic())
    u0 = DynamicsVariables(sim, randn(1,1), randn(1,1))
    run_dynamics(sim, (0.0, 10.0), u0; dt=0.1, reduction, output)

    h5open("test.h5") do fid
        @test Array(fid["trajectory_1"]["OutputPosition"]) isa Array{<:Real,3}
        @test Array(fid["trajectory_1"]["OutputTotalEnergy"]) isa Vector{<:Real}
    end
    rm("test.h5")
end
