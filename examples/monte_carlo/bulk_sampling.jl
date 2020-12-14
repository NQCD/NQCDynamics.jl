using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.IO
using PyCall
using Unitful
using Plots

build = pyimport("ase.build")
ase = pyimport("ase")

pd = build.bulk("Pd", "fcc", a=3.89, cubic=true)
pd += ase.Atoms("H", positions=[[1.5, 1.5, 1.5]])

atoms, positions = extract_parameters_and_positions(pd)
model = Models.PdH(atoms.atom_types, atoms.cell, 10.0)

Δ = Dict([(:Pd, 0.5), (:H, 1.0)])
sys = System{MonteCarlo}(atoms, model, 3000u"K", Δ; passes=2000, fix=[1])

output = InitialConditions.run_monte_carlo_sampling(sys, positions)

@show output.acceptance
write_trajectory("sampling.xyz", output.R, atoms)
plot(output.energy)