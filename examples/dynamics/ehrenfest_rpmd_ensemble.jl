using NonadiabaticMolecularDynamics
using Plots
using Unitful

n_beads = 5
atoms = Atoms([:H])

sim1 = Simulation{Ehrenfest}(atoms, TullyModelOne(); DoFs=1)
sim2 = Simulation{Ehrenfest}(atoms, TullyModelTwo(); DoFs=1)
sim3 = RingPolymerSimulation{Ehrenfest}(atoms, TullyModelOne(), n_beads; DoFs=1) #, temperature=10u"K"
sim4 = RingPolymerSimulation{Ehrenfest}(atoms, TullyModelTwo(), n_beads; DoFs=1)

r = fill(-5.0, sim1.DoFs, length(sim1.atoms))
rr = RingPolymerArray(fill(-5.0, sim3.DoFs, length(sim3.atoms), n_beads))

v1 = fill(8.9, sim1.DoFs, length(sim1.atoms)) ./ sim1.atoms.masses[1]
v2 = fill(16, sim2.DoFs, length(sim2.atoms)) ./ sim2.atoms.masses[1]
vv1 = RingPolymerArray(fill(8.9, sim3.DoFs, length(sim3.atoms), n_beads) ./ sim3.atoms.masses[1])
vv2 = RingPolymerArray(fill(16, sim4.DoFs, length(sim4.atoms), n_beads) ./ sim4.atoms.masses[1])

output1 = Ensembles.OutputDiabaticPopulation(sim1)
output2 = Ensembles.OutputDiabaticPopulation(sim2)
output3 = Ensembles.OutputDiabaticPopulation(sim3)
output4 = Ensembles.OutputDiabaticPopulation(sim4)

distribution1 = InitialConditions.DynamicalDistribution(v1,r,(1,1);state=1,type=:diabatic)
distribution2 = InitialConditions.DynamicalDistribution(v2,r,(1,1);state=1,type=:diabatic)
distribution3 = InitialConditions.DynamicalDistribution(vv1,rr,(1,1);state=1,type=:diabatic)
distribution4 = InitialConditions.DynamicalDistribution(vv2,rr,(1,1);state=1,type=:diabatic)

selection1 = Ensembles.RandomSelection(distribution1)
selection2 = Ensembles.RandomSelection(distribution2)
selection3 = Ensembles.RandomSelection(distribution3)
selection4 = Ensembles.RandomSelection(distribution4)

reduction = Ensembles.MeanReduction()

rng=0:10:2500
  
#model 1
#solution1 = Ensembles.run_ensemble(sim1, (0.0, 2500.0), selection1; trajectories=1e3,
#    output=output1, reduction=reduction, saveat=10.0)
#plot1 = plot(rng, [p[1] for p in solution1.u], title="Tully model 1", label="Ehrenfest P1", legend=:right)
#plot!(rng, [p[2] for p in solution1.u], label="Ehrenfest P2")

solution3 = Ensembles.run_ensemble(sim3, (0.0, 2500.0), selection3; trajectories=1e3,
    output=output3, reduction=reduction, saveat=10.0)
plot!(rng, [p[1] for p in solution3.u], label="Ehrenfest RPMD P1")
plot!(rng, [p[2] for p in solution3.u], label="Ehrenfest RPMD P2")

#model 2
solution2 = Ensembles.run_ensemble(sim2, (0.0, 2500.0), selection2; trajectories=1e3,
    output=output2, reduction=reduction, saveat=10.0)
plot2 = plot(rng, [p[1] for p in solution2.u], title="Tully model 2", label="Ehrenfest P1", legend=:right)
plot!(rng, [p[2] for p in solution2.u], label="Ehrenfest P2")

solution4 = Ensembles.run_ensemble(sim4, (0.0, 2500.0), selection4; trajectories=1e3,
    output=output4, reduction=reduction, saveat=10.0)
plot!(rng, [p[1] for p in solution4.u], label="Ehrenfest RPMD P1")
plot!(rng, [p[2] for p in solution4.u], label="Ehrenfest RPMD P2")

plot(plot1, plot2, layout = (1, 2), size=(1000,400))
