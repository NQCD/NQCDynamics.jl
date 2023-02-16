using Test
using NQCDynamics
using OrdinaryDiffEq: Tsit5
using Random: seed!

@testset "Ehrenfest" begin
    sim = RingPolymerSimulation{Ehrenfest}(Atoms(2000), TullyModelOne(), 10; temperature=1e-3)
    u = DynamicsVariables(sim, fill(10/2000, size(sim)), fill(-5, size(sim)) .+ randn(size(sim)), PureState(1))
    dt = 0.1

    sol = run_dynamics(sim, (0, 2000.0), u; output=OutputDynamicsVariables, algorithm=DynamicsMethods.IntegrationAlgorithms.BCBwithTsit5(Tsit5()), saveat=0:10:2000, dt=dt)
    sol1 = run_dynamics(sim, (0, 2000.0), u; output=OutputDynamicsVariables, algorithm=Tsit5(), saveat=0:10:2000, abstol=1e-6, reltol=1e-6)

    @test sol[:OutputDynamicsVariables] ≈ sol1[:OutputDynamicsVariables] rtol=1e-3
end

@testset "FSSH" begin
    sim = RingPolymerSimulation{FSSH}(Atoms(2000), TullyModelOne(), 10; temperature=1e-3)
    u = DynamicsVariables(sim, fill(20/2000, size(sim)), fill(-5, size(sim)) .+ randn(size(sim)), PureState(1))
    dt = 0.1

    seed!(1)
    sol = run_dynamics(sim, (0, 2000.0), u; output=OutputDynamicsVariables, algorithm=DynamicsMethods.IntegrationAlgorithms.BCBwithTsit5(Tsit5()), saveat=0:10:2000, dt=dt)
    seed!(1)
    sol1 = run_dynamics(sim, (0, 2000.0), u; output=OutputDynamicsVariables, algorithm=Tsit5(), saveat=0:10:2000, dt=dt, adaptive=false)

    @test sol[:OutputDynamicsVariables] ≈ sol1[:OutputDynamicsVariables] rtol=1e-3
end
