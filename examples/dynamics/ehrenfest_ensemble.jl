using NonadiabaticMolecularDynamics
using Plots
using Unitful
using UnitfulAtomic
using Distributions

model = TullyModelOne()
atoms = Atoms([:H])
k = 8.9 
sim = Simulation{Ehrenfest}(atoms, model; DoFs=1)
output1 = Ensembles.OutputDiabaticPopulation(sim)
#output2 = Ensembles.OutputAdiabaticPopulation(sim)
v = k/sim.atoms.masses[1]
r = Normal(-10)
distribution = InitialConditions.DynamicalDistribution(v, r, (1,1); state=2, type=:diabatic)
selection = Ensembles.RandomSelection(distribution)
reduction = Ensembles.MeanReduction()

solution = Ensembles.run_ensemble(sim, (0.0, 3000.0), selection; trajectories=1e3,
    output=output1, reduction=reduction, saveat=10.0)

plot(0:10:3000, [p[1] for p in solution.u], legend=false)
plot!(0:10:3000, [p[2] for p in solution.u])
ylabel!("Population difference")
xlabel!("Time")


