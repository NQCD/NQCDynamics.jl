using NonadiabaticMolecularDynamics
using Plots
using Unitful
using Distributions

n_beads = 5
atoms = Atoms(2000)
k1=8.9 #au 0.0198025 0.064
k2=16
e1=(k1^2)/(2*atoms.masses[1])
e2=(k2^2)/(2*atoms.masses[1])

#e1 = k1 #0.0001u"K"
#e2 = k2 #0.0001u"K"

width = 1/sqrt(0.5)

sim1 = Simulation{FSSH}(atoms, TullyModelOne(); DoFs=1)
sim2 = Simulation{FSSH}(atoms, TullyModelTwo(); DoFs=1)
sim3 = RingPolymerSimulation{FSSH}(atoms, TullyModelOne(), n_beads; DoFs=1, temperature=e1)
sim4 = RingPolymerSimulation{FSSH}(atoms, TullyModelTwo(), n_beads; DoFs=1, temperature=e2)

traj_num = 20

r1 = Normal(-5.0, width)
#r1 = collect(range(-10, -5, length=20))
#r1 = Normal(-4.0) #fill(Normal(-5.0), sim1.DoFs, length(sim1.atoms))
r2 = Normal(-10.0, width)
#rr =  Normal(-4.0) #RingPolymerArray(fill(Normal(-5.0), sim3.DoFs, length(sim3.atoms), n_beads))

v1 = fill(8.9, sim1.DoFs, length(sim1.atoms)) ./ sim1.atoms.masses[1]
v2 = fill(16, sim2.DoFs, length(sim2.atoms)) ./ sim2.atoms.masses[1]
#vv1 = fill(8.9, sim3.DoFs, length(sim3.atoms), n_beads) ./ sim3.atoms.masses[1]
#vv2 = fill(16, sim4.DoFs, length(sim4.atoms), n_beads) ./ sim4.atoms.masses[1]
#x = rand(Normal(0, 0.4), 100)
#rand(Truncated(Normal(0, 0.4), sim3.DoFs, length(sim3.atoms)), n_beads)
vv1 = fill(8.9, sim3.DoFs, length(sim3.atoms), n_beads) ./ sim3.atoms.masses[1]
vv2 = fill(16, sim4.DoFs, length(sim4.atoms), n_beads) ./ sim4.atoms.masses[1]
#vv1 =  (fill(8.9, sim3.DoFs, length(sim3.atoms), n_beads)+ randn(sim3.DoFs, length(sim3.atoms), n_beads)) ./ sim3.atoms.masses[1]
#vv2 =  (fill(16, sim4.DoFs, length(sim4.atoms), n_beads)+ randn(sim4.DoFs, length(sim4.atoms), n_beads)) ./ sim4.atoms.masses[1]

output1 = Ensembles.OutputDiabaticPopulation(sim1)
output2 = Ensembles.OutputDiabaticPopulation(sim2)
output3 = Ensembles.OutputDiabaticPopulation(sim3)
output4 = Ensembles.OutputDiabaticPopulation(sim4)

distribution1 = InitialConditions.DynamicalDistribution(v1,r1,(1,1);state=1,type=:diabatic)
distribution2 = InitialConditions.DynamicalDistribution(v2,r2,(1,1);state=1,type=:diabatic)
distribution3 = InitialConditions.DynamicalDistribution(vv1,r1,(1,1, n_beads);state=1,type=:diabatic)
distribution4 = InitialConditions.DynamicalDistribution(vv2,r2,(1,1, n_beads);state=1,type=:diabatic)

selection1 = Ensembles.OrderedSelection(distribution1)
selection2 = Ensembles.OrderedSelection(distribution2)
selection3 = Ensembles.OrderedSelection(distribution3)
selection4 = Ensembles.OrderedSelection(distribution4)

#reduction = Ensembles.MeanReduction()
reduction = Ensembles.SumReduction([zeros(2) for _ in 1:251])

rng=0:10:3000
rng2=1:301
#model 1
@time solution1 = Ensembles.run_ensemble_standard_output(sim1, (0.0, 3000.0), selection1; trajectories=traj_num,
    output=(:population, :state), saveat=rng)
plt = plot(legend=false)
global avg = [zeros(2) for _ in rng] # zero(solution1[1].population)
for s in solution1
    #global avg += s.population
    plot!(s, :population, label="FSSH", linestyle=:dash)
end
#avg = avg/traj_num
#plot!(rng, [p[1] for p in avg], title="Tully model 1", label="FSSH P1", linestyle=:dash)
#plot!(rng, [p[2] for p in avg], title="Tully model 1", label="FSSH P1", linestyle=:dash)

@time solution3 = Ensembles.run_ensemble_standard_output(sim3, (0.0, 3000.0), selection3; trajectories=traj_num,
    output=(:population, :state), saveat=rng)
global avg = [zeros(2) for _ in rng] #zero(solution3[1].population)
for s in solution3
    #plot!(rng, [p[1]-p[2] for p in s.population], title="Tully model 1", label="FSSH P1")
    #global avg += s.population
    plot!(s, :population, label="FSSH P2")
    #plot!(rng, [p[2] for p in s.population], label="FSSH P2")
end
#avg = avg/traj_num
#plot!(rng, [p[1] for p in avg], title="Tully model 1", label="FSSH P1")
#plot!(rng, [p[2] for p in avg], title="Tully model 1", label="FSSH P1")

plt
# #model 2
# @time solution2 = Ensembles.run_ensemble_standard_output(sim2, (0.0, 2500.0), selection2; trajectories=1e3,
#     output=(:population), saveat=10.0)
# plot2 = plot(rng, [p[1] for p in solution2.u], title="Tully model 2", label="FSSH P1", legend=:right)
# plot!(rng, [p[2] for p in solution2.u], label="FSSH P2")

# @time solution4 = Ensembles.run_ensemble_standard_output(sim4, (0.0, 2500.0), selection4; trajectories=1e3,
#     output=(:population), saveat=10.0)
# plot!(rng, [p[1] for p in solution4.u], label="FSSH RPMD P1")
# plot!(rng, [p[2] for p in solution4.u], label="FSSH RPMD P2")

# plot(plot1, plot2, layout = (1, 2), size=(1000,400))
