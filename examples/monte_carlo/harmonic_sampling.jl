using NonadiabaticMolecularDynamics
using Unitful
using Plots

atoms = Atoms.AtomicParameters(Atoms.InfiniteCell(), fill(:H, 5))
model = Models.Harmonic(1, 1, 0.1)

Δ = Dict([(:H, 0.3)])
sys = System{MonteCarlo}(atoms, model, 300u"K", Δ, 1; passes=2000)

output = InitialConditions.run_monte_carlo_sampling(sys, zeros(1, n_atoms(sys)))

@show output.acceptance
plot(vcat(output.R...))