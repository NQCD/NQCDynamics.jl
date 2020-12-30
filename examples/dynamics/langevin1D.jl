using NonadiabaticMolecularDynamics
using Plots
using Unitful

atoms = Atoms{Float64}([:H, :H])
cell = InfiniteCell{Float64}()
model = Models.Harmonic()
sim = Simulation(1, 300u"K", cell, atoms, model, Dynamics.Langevin(1))

z = Dynamics.Phasespace(zeros(1, length(atoms)), zeros(1, length(atoms)))

solution = Dynamics.run_trajectory(z, (0.0, 50.0), sim)
plot(solution, vars=range(atoms))