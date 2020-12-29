using NonadiabaticMolecularDynamics
using Plots
using Random

atoms = Atoms{Float64}(vcat(fill(:H, 5), fill(:C, 2)))
cell = InfiniteCell{Float64}()
model = Models.Harmonic()
sim = Simulation(1, 1, cell, atoms, model, Dynamics.Classical())

z = Phasespace(randn(1, length(atoms)), randn(1, length(atoms)))

solution = Dynamics.run_trajectory(z, (0.0, 1e3), sim)
plot(solution, vars=range(atoms))
