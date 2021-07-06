using NonadiabaticMolecularDynamics
using Plots
using Unitful
n_beads = 4

atoms = Atoms(2000)
sim1 = Simulation{FSSH}(atoms, TullyModelOne(); DoFs=1)
sim2 = Simulation{FSSH}(atoms, TullyModelTwo(); DoFs=1)
k1=8.9
k2=16
r = fill(-5.0, sim1.DoFs, length(sim1.atoms))
v1 = fill(k1, sim1.DoFs, length(sim1.atoms)) ./ sim1.atoms.masses[1]
v2 = fill(k2, sim2.DoFs, length(sim2.atoms)) ./ sim2.atoms.masses[1]
z1 = SurfaceHoppingVariables(v1, r, 2, 1)
z2 = SurfaceHoppingVariables(v2, r, 2, 1)

e1=(k1^2)/(2*atoms.masses[1])
e2=(k2^2)/(2*atoms.masses[1])

sim3 = RingPolymerSimulation{FSSH}(atoms, TullyModelOne(), n_beads; DoFs=1, temperature=e1)
sim4 = RingPolymerSimulation{FSSH}(atoms, TullyModelTwo(), n_beads; DoFs=1, temperature=e2)
rr = RingPolymerArray(fill(-5.0, sim3.DoFs, length(sim3.atoms), n_beads))
vv1 = RingPolymerArray(fill(8.9, sim3.DoFs, length(sim3.atoms), n_beads) ./ sim3.atoms.masses[1])
vv2 = RingPolymerArray(fill(16, sim4.DoFs, length(sim4.atoms), n_beads) ./ sim4.atoms.masses[1])
z3 = SurfaceHoppingVariables(vv1, rr, 2, 1)
z4 = SurfaceHoppingVariables(vv2, rr, 2, 1)

#model 1
#@time solution1 = Dynamics.run_trajectory(z1, (0.0, 2500.0), sim1; output=(:state))
#plot1 = plot(solution1.t, solution1.state, title="Tully model 1", label="FSSH P1", legend=:right)
#plot!(solution1.t, solution1.state, label="FSSH P2")
#plot!(solution1.t, [p[2] for p in solution1.state], label="FSSH P2")

@time solution3 = Dynamics.run_trajectory(z3, (0.0, 2500.0), sim3; output=(:population))

# plot!(solution3.t, [p[1] for p in solution3.state], label="FSSH RPMD P1")
# plot!(solution3.t, [p[2] for p in solution3.state], label="FSSH RPMD P2")
#plot!(solution3.t, solution3.state, label="FSSH P2")
plot1 = plot(solution3, :population)
# #model 2
# @time solution2 = Dynamics.run_trajectory(z2, (0.0, 2500.0), sim2; output=(:population))
# plot2 = plot(solution2.t, [p[1] for p in solution2.population], title="Tully model 2", label="FSSH P1", legend=:right)
# plot!(solution2.t, [p[2] for p in solution2.population], label="FSSH P2")

# @time solution4 = Dynamics.run_trajectory(z4, (0.0, 2500.0), sim4; output=(:population))
# plot!(solution4.t, [p[1] for p in solution4.population], label="FSSH RPMD P1")
# plot!(solution4.t, [p[2] for p in solution4.population], label="FSSH RPMD P2")

# plot(plot1, plot2, layout = (1, 2), size=(1000,400))
