using NonadiabaticMolecularDynamics
using PyCall
using Unitful
using Plots
import JuLIP
using ASE

build = pyimport("ase.build")
ase = pyimport("ase")

pd = build.bulk("Pd", "fcc", a=3.89, cubic=true)
pd += ase.Atoms("H", positions=[[1.5, 1.5, 1.5]])

cell, atoms, positions = extract_parameters_and_positions(pd)

model = Models.JuLIPModel(atoms, cell, JuLIP.EAM("PdH_hijazi.eam.alloy"))

Δ = Dict([(:Pd, 0.5), (:H, 1.0)])
sim = Simulation(atoms, model; temperature=300u"K", cell=cell)

output = InitialConditions.run_monte_carlo_sampling(sim, positions, Δ, 20, Int[])

@show output.acceptance
write_trajectory("sampling.xyz", cell, atoms, output.R)
plot(output.energy)