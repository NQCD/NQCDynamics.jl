using NonadiabaticMolecularDynamics
using Unitful
using Plots

atoms = Atoms.AtomicParameters(Atoms.InfiniteCell(), fill(:H, 5))
model = Models.Harmonic(1, 1, 0.1)

sys = System(atoms, model, 1; temperature=300u"K")

step_sizes = Dict([(:H, 0.3)])
output = InitialConditions.run_monte_carlo_sampling(sys, zeros(1, n_atoms(sys)), step_sizes; passes=2000)

@show output.acceptance
plot(vcat(output.R...))