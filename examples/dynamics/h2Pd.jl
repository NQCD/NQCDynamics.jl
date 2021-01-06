using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.IO
using PyCall

build = pyimport("ase.build")
ase = pyimport("ase")

slab = build.fcc111("Pd", size=(3,2,3), orthogonal=true)
build.add_adsorbate(slab, "H", 1.5, "ontop")
build.add_adsorbate(slab, "H", 1.5, "bridge")
slab.center(vacuum=10.0, axis=2)

cell, atoms, positions = extract_parameters_and_positions(slab)
model = Models.PdH(atoms.types, cell, 10.0)

sim = Simulation(atoms, model, Dynamics.Classical(); cell=cell)

z = Phasespace(positions, zero(positions))

@time solution = Dynamics.run_trajectory(z, (0.0, 1000.0), sim)
write_trajectory("traj.xyz", cell, atoms, get_positions.(solution.u))
