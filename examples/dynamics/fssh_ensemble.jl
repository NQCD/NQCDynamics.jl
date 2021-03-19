using NonadiabaticMolecularDynamics
using Unitful
using UnitfulAtomic
using Distributions

model = Models.TullyModelOne()

atoms = NonadiabaticMolecularDynamics.Atoms([:H])

k = 10
sim = Simulation{FSSH}(atoms, model; DoFs=1)

distribution = InitialConditions.DynamicalDistribution(k / atoms.masses[1], Normal(1), (1,1))

@time sol = Dynamics.run_ensemble(distribution, (0.0, 3000.0), sim; trajectories=10, output=(:velocity))

plt = plot()
for i=1:10
    plot!(sol[i], :velocity)
end
plt
