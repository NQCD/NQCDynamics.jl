using NonadiabaticMolecularDynamics
using Unitful
using Plots

atoms = Atoms(fill(:H, 5))
model = Models.Harmonic()

Δ = Dict([(:H, 0.3)])
monte_carlo = InitialConditions.MonteCarlo{Float64}(Δ, length(atoms), 1000, Int[])
sim = Simulation{Classical}(atoms, model; DoFs=1, temperature=300u"K")

output = InitialConditions.run_monte_carlo_sampling(sim, zeros(1, length(atoms)), Δ, 1000)

@show output.acceptance
plot(vcat(output.R...))