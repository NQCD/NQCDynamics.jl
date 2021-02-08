using Test
using NonadiabaticMolecularDynamics

cb, vals = Dynamics.create_energy_saving_callback()

atoms = Atoms(vcat(fill(:H, 5), fill(:C, 2)))
model = Models.Harmonic()
sim = Simulation(atoms, model, Classical())

z = ClassicalDynamicals(randn(3, length(atoms)), randn(3, length(atoms)))

solution = Dynamics.run_trajectory(z, (0.0, 1e3), sim; callback=cb, dt=0.1)

@test solution.t ≈ vals.t
@test Models.energy.(model, get_positions.(solution.u)) ≈ vals.saveval
