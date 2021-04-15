using NonadiabaticMolecularDynamics
using PyCall
using ASE
using Unitful
using Plots
import JuLIP

build = pyimport("ase.build")
ase = pyimport("ase")

slab = build.fcc111("Pd", size=(3,2,3), orthogonal=true)
build.add_adsorbate(slab, "H", 1.5, "ontop")
build.add_adsorbate(slab, "H", 1.5, "bridge")
slab.center(vacuum=10.0, axis=2)

cell, atoms, positions = extract_parameters_and_positions(slab)
model = Models.JuLIPModel(atoms, cell, JuLIP.Potentials.EAM("PdH_hijazi.eam.alloy"))

Δ = Dict([(:Pd, 0.5), (:H, 1.0)])
sim = RingPolymerSimulation{Classical}(atoms, model, 10; quantum_nuclei=[:H], temperature=100u"K", cell=cell)

R = cat([positions for i=1:10]..., dims=3)
output = InitialConditions.run_monte_carlo_sampling(sim, R, Δ, 100; fix=collect(1:6))

@show output.acceptance
write_trajectory("sampling.xyz", cell, atoms, output.R)
plot(output.energy)