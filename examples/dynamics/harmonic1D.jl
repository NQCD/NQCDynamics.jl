using NonadiabaticMolecularDynamics
using Plots
using Random

atoms = Atoms(vcat(fill(:H, 5), fill(:C, 2)))
model = Harmonic()
sim = Simulation{Classical}(atoms, model; DoFs=1)

z = ClassicalDynamicals(zeros(1, length(atoms)), randn(1, length(atoms)))

solution = Dynamics.run_trajectory(z, (0.0, 1e3), sim, dt=0.1, output=(:hamiltonian, :potential))
plot(solution, :hamiltonian, label="hamiltonian", legend=:right)
plot!(solution, :potential, label="potential energy", legend=:right)
