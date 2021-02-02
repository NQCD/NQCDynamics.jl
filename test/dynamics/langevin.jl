using Test
using NonadiabaticMolecularDynamics
using Unitful

@test Dynamics.Langevin{Float64}(1.0) isa Dynamics.Langevin

atoms = Atoms([:H, :H])
model = Models.Harmonic()

sim = Simulation{Langevin}(atoms, model; DoFs=2, temperature=300u"K")

z = Dynamics.Phasespace(zeros(2, length(atoms)), zeros(2, length(atoms)))

solution = Dynamics.run_trajectory(z, (0.0, 50.0), sim)
