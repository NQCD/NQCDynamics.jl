using NonadiabaticMolecularDynamics
using Plots
using Unitful

atoms = Atoms([:H])
sim1 = Simulation{Ehrenfest}(atoms, Models.TullyModelOne(); DoFs=1)
sim2 = Simulation{Ehrenfest}(atoms, Models.TullyModelTwo(); DoFs=1)
r = fill(-5.0, sim1.DoFs, length(sim1.atoms))
v1 = fill(8.9, sim1.DoFs, length(sim1.atoms)) ./ sim1.atoms.masses[1]
v2 = fill(16, sim1.DoFs, length(sim1.atoms)) ./ sim1.atoms.masses[1]
z1 = EhrenfestVariables(v1, r, 2)#, 1)
z2 = EhrenfestVariables(v2, r, 2)#, 1)
solution1 = Dynamics.run_trajectory(z1, (0.0, 2500.0), sim1; output=(:population))
plot1 = plot(solution1.t, [p[1] for p in solution1.population], title="Tully model 1", label="Ehrenfest P1")
plot!(solution1.t, [p[2] for p in solution1.population], label="Ehrenfest P2")

solution2 = Dynamics.run_trajectory(z2, (0.0, 2500.0), sim2; output=(:population))
plot2 = plot(solution2.t, [p[1] for p in solution2.population], title="Tully model 2", label="Ehrenfest P1")
plot!(solution2.t, [p[2] for p in solution2.population], label="Ehrenfest P2")

plot(plot1, plot2, layout = (1, 2), size=(1000,400))