using NonadiabaticMolecularDynamics
using Unitful
using Plots

# atoms = Atoms.AtomicParameters(Atoms.InfiniteCell(), fill(:H, 5))
atoms = Atoms{Float64}(fill(:H, 5))
cell = InfiniteCell{Float64}()
calc = Calculators.Harmonic()
# model = Models.Harmonic(1, 1, 0.1)

Δ = Dict([(:H, 0.3)])
monte_carlo = InitialConditions.MonteCarlo{Float64}(Δ, length(atoms), 1000, Int[])
sim = Simulation(1, 300u"K", cell, atoms, calc, monte_carlo)

output = InitialConditions.run_monte_carlo_sampling(sim, zeros(1, length(sim.atoms)))

@show output.acceptance
plot(vcat(output.R...))