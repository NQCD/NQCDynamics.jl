using NonadiabaticMolecularDynamics
using Unitful
using UnitfulAtomic
using Distributions

model = Models.TullyModelOne()

atoms = NonadiabaticMolecularDynamics.Atoms([:H])

k = 10
sim = Simulation{FSSH}(atoms, model; DoFs=1)

distribution = InitialConditions.PhasespaceDistribution(Normal(1), k, (1,1))

@time sol = Dynamics.run_ensemble(distribution, (0.0, 3000.0), sim; trajectories=10)
