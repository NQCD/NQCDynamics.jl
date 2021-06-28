using Plots: ensure_gradient!
using NonadiabaticMolecularDynamics
using Plots
using Unitful
using UnitfulAtomic
using Distributions

model = TullyModelOne()
model2 = TullyModelTwo()
atoms = Atoms(2000)
k = 8.9 
k2 = 16
sim = Simulation{Ehrenfest}(atoms, model; DoFs=1)
sim2 = Simulation{Ehrenfest}(atoms, model2; DoFs=1)
output1 = Ensembles.OutputDiabaticPopulation(sim)
#output2 = Ensembles.OutputAdiabaticPopulation(sim)
output2 = Ensembles.OutputDiabaticPopulation(sim2)
v = k/sim.atoms.masses[1]
v2 = k2/sim2.atoms.masses[1]
r = Normal(-10)
r2 = Normal(-10)
distribution = InitialConditions.DynamicalDistribution(v, r, (1,1); state=2, type=:diabatic)
selection = Ensembles.RandomSelection(distribution)
reduction = Ensembles.MeanReduction()
distribution2 = InitialConditions.DynamicalDistribution(v2, r2, (1,1); state=2, type=:diabatic)
selection2 = Ensembles.RandomSelection(distribution2)

solution = Ensembles.run_ensemble(sim, (0.0, 3000.0), selection; trajectories=1e3,
    output=output1, reduction=reduction, saveat=10.0)

plot1 = plot(0:10:3000, [p[1] for p in solution.u], title="Tully model 1", legend=false)
plot!(0:10:3000, [p[2] for p in solution.u])


solution2 = Ensembles.run_ensemble(sim2, (0.0, 3000.0), selection2; trajectories=1e3,
    output=output2, reduction=reduction, saveat=10.0)

plot2 = plot(0:10:3000, [p[1] for p in solution2.u], title="Tully model 2", legend=false)
plot!(0:10:3000, [p[2] for p in solution2.u])

plot(plot1, plot2, layout = (1, 2), size=(1000,400))

ylabel!("Population difference")
xlabel!("Time")

