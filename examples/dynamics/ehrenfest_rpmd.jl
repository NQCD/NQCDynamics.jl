using NonadiabaticMolecularDynamics
using Plots
using Unitful
n_beads = 5

atoms = Atoms([:H])
sim1 = Simulation{Ehrenfest}(atoms, TullyModelOne(); DoFs=1)
sim2 = Simulation{Ehrenfest}(atoms, TullyModelTwo(); DoFs=1)
r = fill(-5.0, sim1.DoFs, length(sim1.atoms))
v1 = fill(8.9, sim1.DoFs, length(sim1.atoms)) ./ sim1.atoms.masses[1]
v2 = fill(16, sim2.DoFs, length(sim2.atoms)) ./ sim2.atoms.masses[1]
z1 = EhrenfestVariables(v1, r, 2, 1)
z2 = EhrenfestVariables(v2, r, 2, 1)

sim3 = RingPolymerSimulation{Ehrenfest}(atoms, TullyModelOne(), n_beads; DoFs=1) #, temperature=10u"K"
sim4 = RingPolymerSimulation{Ehrenfest}(atoms, TullyModelTwo(), n_beads; DoFs=1)
rr = RingPolymerArray(fill(-5.0, sim3.DoFs, length(sim3.atoms), n_beads))
vv1 = RingPolymerArray(fill(8.9, sim3.DoFs, length(sim3.atoms), n_beads) ./ sim3.atoms.masses[1])
vv2 = RingPolymerArray(fill(16, sim4.DoFs, length(sim4.atoms), n_beads) ./ sim4.atoms.masses[1])
z3 = EhrenfestVariables(vv1, rr, 2, 1)
z4 = EhrenfestVariables(vv2, rr, 2, 1)

#model 1
solution1 = Dynamics.run_trajectory(z1, (0.0, 2500.0), sim1; output=(:population))
plot1 = plot(solution1.t, [p[1] for p in solution1.population], title="Tully model 1", label="Ehrenfest P1", legend=:right)
plot!(solution1.t, [p[2] for p in solution1.population], label="Ehrenfest P2")

solution3 = Dynamics.run_trajectory(z3, (0.0, 2500.0), sim3; output=(:population))
plot!(solution3.t, [p[1] for p in solution3.population], label="Ehrenfest RPMD P1")
plot!(solution3.t, [p[2] for p in solution3.population], label="Ehrenfest RPMD P2")

#model 2
solution2 = Dynamics.run_trajectory(z2, (0.0, 2500.0), sim2; output=(:population))
plot2 = plot(solution2.t, [p[1] for p in solution2.population], title="Tully model 2", label="Ehrenfest P1", legend=:right)
plot!(solution2.t, [p[2] for p in solution2.population], label="Ehrenfest P2")

solution4 = Dynamics.run_trajectory(z4, (0.0, 2500.0), sim4; output=(:population))
plot!(solution4.t, [p[1] for p in solution4.population], label="Ehrenfest RPMD P1")
plot!(solution4.t, [p[2] for p in solution4.population], label="Ehrenfest RPMD P2")

plot(plot1, plot2, layout = (1, 2), size=(1000,400))
