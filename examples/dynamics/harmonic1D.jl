using NonadiabaticMolecularDynamics
using Plots
using Random

atoms = Atoms(vcat(fill(:H, 5), fill(:C, 2)))
model = Models.Harmonic()
sim = Simulation{Classical}(atoms, model; DoFs=1)

z = ClassicalDynamicals(zeros(1, length(atoms)), randn(1, length(atoms)))

solution = Dynamics.run_trajectory(z, (0.0, 1e3), sim, dt=0.1, output=(:energy, :potential))
plot(solution, :energy)
plot!(solution, :potential)
