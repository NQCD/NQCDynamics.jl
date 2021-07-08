using NonadiabaticMolecularDynamics
using Plots
using Unitful
using Distributions

n_beads = 4

atoms = Atoms(2000)
sim1 = Simulation{Ehrenfest}(atoms, TullyModelOne(); DoFs=1)
sim2 = Simulation{Ehrenfest}(atoms, TullyModelTwo(); DoFs=1)
k1=8.9
k2=16
r = fill(-5.0, sim1.DoFs, length(sim1.atoms))
v1 = fill(k1, sim1.DoFs, length(sim1.atoms)) ./ sim1.atoms.masses[1]
v2 = fill(k2, sim2.DoFs, length(sim2.atoms)) ./ sim2.atoms.masses[1]
z1 = EhrenfestVariables(v1, r, 2, 1)
z2 = EhrenfestVariables(v2, r, 2, 1)

e1=(k1^2)/(2*atoms.masses[1])
e2=(k2^2)/(2*atoms.masses[1])

#e1 = k1
#e2 = k2
normal = Normal(-5.0)

sim3 = RingPolymerSimulation{Ehrenfest}(atoms, TullyModelOne(), n_beads; DoFs=1, temperature=e1)
sim4 = RingPolymerSimulation{Ehrenfest}(atoms, TullyModelTwo(), n_beads; DoFs=1, temperature=e2)
rr = RingPolymerArray(rand(normal, sim3.DoFs, length(sim3.atoms), n_beads))
vv1 =  RingPolymerArray((fill(8.9, sim3.DoFs, length(sim3.atoms), n_beads)) ./ sim3.atoms.masses[1])
vv2 =  RingPolymerArray((fill(16, sim4.DoFs, length(sim4.atoms), n_beads)) ./ sim4.atoms.masses[1])
z3 = EhrenfestVariables(vv1, rr, 2, 1)
z4 = EhrenfestVariables(vv2, rr, 2, 1)

#model 1
@time solution1 = Dynamics.run_trajectory(z1, (0.0, 2500.0), sim1; output=(:population, :quantum_subsystem),saveat=0:1:2500, reltol=1e-5, abstol=1e-8)
plot1 = plot(solution1.t, [p[1] for p in solution1.population], title="Tully model 1", label="Ehrenfest P1", legend=:right)
plot!(solution1.t, [p[2] for p in solution1.population], label="Ehrenfest P2")

#plot1 = plot(solution1.t, [p[1] for p in solution1.population], title="Tully model 1", label="Ehrenfest P1", legend=:right)
#plot(solution1, :quantum_subsystem)

@time solution3 = Dynamics.run_trajectory(z3, (0.0, 2500.0), sim3; output=(:population, :quantum_subsystem),saveat=0:1:2500, reltol=1e-5, abstol=1e-8)
plot!(solution3.t, [p[1] for p in solution3.population], label="Ehrenfest RPMD P1")
plot!(solution3.t, [p[2] for p in solution3.population], label="Ehrenfest RPMD P2")









#plot!(solution3.t, [p[1] for p in solution3.position], label="Ehrenfest RPMD P1")

#plot!(solution3, :quantum_subsystem)

#model 2
# @time solution2 = Dynamics.run_trajectory(z2, (0.0, 2500.0), sim2; output=(:population, :position),saveat=0:1:2500, reltol=1e-5, abstol=1e-8)
# plot2 = plot(solution2.t, [p[1] for p in solution2.population], title="Tully model 2", label="Ehrenfest P1", legend=:right)
# plot!(solution2.t, [p[2] for p in solution2.population], label="Ehrenfest P2")

# @time solution4 = Dynamics.run_trajectory(z4, (0.0, 2500.0), sim4; output=(:population, :position),saveat=0:1:2500, reltol=1e-5, abstol=1e-8)
# plot!(solution4.t, [p[1] for p in solution4.population], label="Ehrenfest RPMD P1")
# plot!(solution4.t, [p[2] for p in solution4.population], label="Ehrenfest RPMD P2")

# plot(plot1, plot2, layout = (1, 2), size=(1000,400))
