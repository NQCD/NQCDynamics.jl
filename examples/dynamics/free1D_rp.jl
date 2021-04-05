using NonadiabaticMolecularDynamics
using Plots
using Random
using Unitful

atoms = Atoms([:H, :C])
model = Models.Free()
sim = RingPolymerSimulation{Classical}(atoms, model, 5; DoFs=1, quantum_nuclei=[:H], temperature=100u"K")

r = RingPolymerArray(randn(1, 2, 5), quantum=sim.beads.quantum_atoms)
v = RingPolymerArray(zeros(1, 2, 5), quantum=sim.beads.quantum_atoms)
z = ClassicalDynamicals(v, r)

@time solution = Dynamics.run_trajectory(z, (0.0, 1e3), sim; dt=0.1, output=(:velocity))
plot(solution, :velocity)