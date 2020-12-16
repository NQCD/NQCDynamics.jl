using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.IO
using PyCall
using Unitful
using Plots

build = pyimport("ase.build")
ase = pyimport("ase")

slab = build.fcc111("Pd", size=(3,2,3), orthogonal=true)
build.add_adsorbate(slab, "H", 1.5, "ontop")
build.add_adsorbate(slab, "H", 1.5, "bridge")
slab.center(vacuum=10.0, axis=2)

atoms, positions = extract_parameters_and_positions(slab)
model = Models.PdH(atoms.atom_types, atoms.cell, 10.0)

Δ = Dict([(:Pd, 0.5), (:H, 0.5)])
sys = System{MonteCarlo}(atoms, model, 300u"K", Δ; passes=1000, fix=collect(1:6))

output = InitialConditions.run_monte_carlo_sampling(sys, positions)

@show output.acceptance
write_trajectory("sampling.xyz", output.R, atoms)
plot(output.energy)