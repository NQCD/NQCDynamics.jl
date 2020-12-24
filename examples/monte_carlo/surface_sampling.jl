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

cell, atoms, positions = extract_parameters_and_positions(slab)
model = Models.PdH(atoms.types, cell, 10.0)

Δ = Dict([(:Pd, 0.5), (:H, 0.5)])
monte_carlo = InitialConditions.MonteCarlo{Float64}(Δ, length(atoms), 100, collect(1:6))
sim = RingPolymerSimulation(3, 100u"K", cell, atoms, model, monte_carlo, 10, [:H])

R = cat([positions for i=1:10]..., dims=3)
output = InitialConditions.run_monte_carlo_sampling(sim, R)

@show output.acceptance
write_trajectory("sampling.xyz", cell, atoms, output.R)
plot(output.energy)