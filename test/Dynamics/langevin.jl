using Test
using NQCDynamics

atoms = Atoms([:H])
model = Harmonic()

sim = Simulation{Langevin}(atoms, model; temperature=1, γ=1)

z = DynamicsVariables(sim, randn(1,1), randn(1,1))

solution = run_dynamics(sim, (0.0, 500.0), z; output=OutputDynamicsVariables, dt=1)
