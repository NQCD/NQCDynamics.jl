using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.IO
using PyCall
using Unitful
using Plots

build = pyimport("ase.build")
ase = pyimport("ase")

slab = build.fcc111("Pd", size=(3,2,4), orthogonal=true)
build.add_adsorbate(slab, "H", 1.5, "ontop")
build.add_adsorbate(slab, "H", 1.5, "bridge")
slab.center(vacuum=10.0, axis=2)

atoms, positions = extract_parameters_and_positions(slab)
model = Models.PdH(atoms.atom_types, atoms.cell, 10.0)

sys = System(atoms, model; temperature=300u"K")

step_sizes = Dict([(:Pd, 0.5), (:H, 0.5)])
output = InitialConditions.run_monte_carlo_sampling(sys, positions, step_sizes; passes=1000, fix=collect(1:6))

@show output.acceptance
write_trajectory("sampling.xyz", output.R, atoms)
plot(output.energy)