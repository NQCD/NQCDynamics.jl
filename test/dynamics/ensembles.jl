using NonadiabaticMolecularDynamics
using Test

atoms = Atoms{Float64}([:H])
model = Models.Harmonic()
sim = Simulation(atoms, model, Dynamics.Classical(); DoFs=1)

positions = [randn(1, length(atoms)) for i=1:10]
velocities = [randn(1, length(atoms)) for i=1:10]
distribution = InitialConditions.DynamicalDistribution(positions, velocities, (1, 1))

solution = Dynamics.run_ensemble(distribution, (0.0, 1e3), sim; trajectories=100, dt=1)
