using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.IO
using PyCall
using Unitful
using Plots

build = pyimport("ase.build")
ase = pyimport("ase")

pd = build.bulk("Pd", "fcc", a=3.89, cubic=true)
pd += ase.Atoms("H", positions=[[1.5, 1.5, 1.5]])

cell, atoms, positions = extract_parameters_and_positions(pd)
model = Models.PdH(atoms.types, cell, 10.0)

Δ = Dict([(:Pd, 0.5), (:H, 1.0)])
monte_carlo = InitialConditions.MonteCarlo{Float64}(Δ, length(atoms), 2000, [1])
sim = Simulation(atoms, model, Dynamics.Classical(); temperature=300u"K", cell=cell)

output = InitialConditions.run_monte_carlo_sampling(sim, monte_carlo, positions)

@show output.acceptance
write_trajectory("sampling.xyz", cell, atoms, output.R)
plot(output.energy)