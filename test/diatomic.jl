using Test
using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.InitialConditions

atoms = Atoms{Float64}([:N, :O])
model = Models.Harmonic()
sim = Simulation(atoms, model, Dynamics.Classical())

Diatomic.calculate_diatomic_energy(sim, 10.0)