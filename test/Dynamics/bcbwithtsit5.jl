using Test
using NonadiabaticMolecularDynamics
using OrdinaryDiffEq: Tsit5
using Random: seed!

@testset "Ehrenfest" begin
    sim = RingPolymerSimulation{Ehrenfest}(Atoms(2000), TullyModelOne(), 10; temperature=1e-3)
    u = DynamicsVariables(sim, fill(10/2000, size(sim)), fill(-5, size(sim)) .+ randn(size(sim)), SingleState(1))
    dt = 0.1

    sol = run_trajectory(u, (0, 2000.0), sim; algorithm=DynamicsMethods.IntegrationAlgorithms.BCBwithTsit5(), saveat=0:10:2000, dt=dt)
    sol1 = run_trajectory(u, (0, 2000.0), sim; algorithm=Tsit5(), saveat=0:10:2000, abstol=1e-6, reltol=1e-6)

    @test sol.u ≈ sol1.u rtol=1e-3
end

@testset "FSSH" begin
    sim = RingPolymerSimulation{FSSH}(Atoms(2000), TullyModelOne(), 10; temperature=1e-3)
    u = DynamicsVariables(sim, fill(20/2000, size(sim)), fill(-5, size(sim)) .+ randn(size(sim)), SingleState(1))
    dt = 0.1

    seed!(1)
    sol = run_trajectory(u, (0, 2000.0), sim; algorithm=DynamicsMethods.IntegrationAlgorithms.BCBwithTsit5(), saveat=0:10:2000, dt=dt)
    seed!(1)
    sol1 = run_trajectory(u, (0, 2000.0), sim; algorithm=Tsit5(), saveat=0:10:2000, dt=dt, adaptive=false)

    @test sol.u ≈ sol1.u rtol=1e-3
end
