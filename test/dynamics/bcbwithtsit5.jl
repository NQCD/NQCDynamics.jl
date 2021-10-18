using Test
using NonadiabaticMolecularDynamics
using OrdinaryDiffEq: Tsit5

sim = RingPolymerSimulation{Ehrenfest}(Atoms(2000), TullyModelOne(), 10; temperature=1e-3)
u = DynamicsVariables(sim, fill(10/2000, size(sim)), fill(-5, size(sim)) .+ randn(size(sim)), 1)
dt = 1.0

sol = run_trajectory(u, (0, 2000.0), sim; algorithm=DynamicsMethods.IntegrationAlgorithms.BCBwithTsit5(), output=(:population), saveat=0:10:2000, dt=dt)
sol1 = run_trajectory(u, (0, 2000.0), sim; algorithm=Tsit5(), output=(:population), saveat=0:10:2000)

@test sol.population ≈ sol1.population rtol=1e-2
