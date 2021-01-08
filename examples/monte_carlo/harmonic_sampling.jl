using NonadiabaticMolecularDynamics
using Unitful
using Plots

atoms = Atoms{Float64}(fill(:H, 5))
model = Models.Harmonic()

Δ = Dict([(:H, 0.3)])
monte_carlo = InitialConditions.MonteCarlo{Float64}(Δ, length(atoms), 1000, Int[])
sim = Simulation(atoms, model, Dynamics.Classical(); DoFs=1, temperature=300u"K")

output = InitialConditions.run_monte_carlo_sampling(sim, monte_carlo, zeros(1, length(sim.atoms)))

@show output.acceptance
plot(vcat(output.R...))