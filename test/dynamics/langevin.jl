using Test
using NonadiabaticMolecularDynamics

atoms = Atoms([:H])
model = Models.Harmonic()

sim = Simulation{Langevin}(atoms, model; DoFs=1, temperature=1, γ=1)

z = Dynamics.ClassicalDynamicals(randn(sim.DoFs, length(atoms)), randn(sim.DoFs, length(atoms)))

solution = Dynamics.run_trajectory(z, (0.0, 500.0), sim; dt=1)
