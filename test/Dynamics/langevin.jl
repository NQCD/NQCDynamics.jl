using Test
using NonadiabaticMolecularDynamics

atoms = Atoms([:H])
model = Harmonic()

sim = Simulation{Langevin}(atoms, model; temperature=1, γ=1)

z = DynamicsVariables(sim, randn(1,1), randn(1,1))

solution = run_trajectory(z, (0.0, 500.0), sim; dt=1)
