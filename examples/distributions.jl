using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.InitialConditions
using Distributions
using Unitful
using UnitfulAtomic
using Plots

atoms = Atoms([:H])
temperature = austrip(100u"K")
boltzmann = Normal(temperature/atoms.masses[1])

Δ = Dict([(:H, 0.3)])
monte_carlo = InitialConditions.MonteCarlo{Float64}(Δ, length(atoms), 1000, Int[])
model = Models.Harmonic()
sim = Simulation(atoms, model; DoFs=1, temperature=100u"K")

output = InitialConditions.run_monte_carlo_sampling(sim, monte_carlo, zeros(1, length(sim.atoms)))

# @show output.acceptance
# plot(vcat(output.R...))

d = PhasespaceDistribution(output.R, boltzmann, (1,1)) # final argument is size 
display(rand(d))
