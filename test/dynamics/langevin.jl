using Test
using NonadiabaticMolecularDynamics

atoms = Atoms([:H])
model = NonadiabaticModels.Harmonic()

sim = Simulation{Langevin}(atoms, model; DoFs=1, temperature=1, γ=1)

z = ComponentVector(v=randn(sim.DoFs, length(atoms)), r=randn(sim.DoFs, length(atoms)))

solution = Dynamics.run_trajectory(z, (0.0, 500.0), sim; dt=1)
