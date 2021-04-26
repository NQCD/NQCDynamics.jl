using NonadiabaticMolecularDynamics
using Plots
using Unitful

atoms = Atoms([:H])
sim = Simulation{FSSH}(atoms, Models.TullyModelOne(); DoFs=1)

r = fill(-5.0, sim.DoFs, length(sim.atoms))
v = fill(8.9, sim.DoFs, length(sim.atoms)) ./ sim.atoms.masses[1]
z = SurfaceHoppingVariables(v, r, 2, 1)

solution = Dynamics.run_trajectory(z, (0.0, 2500.0), sim; output=(:population, :state))

plot(solution.t, [p[1] for p in solution.population], label="FSSH P1")
plot!(solution.t, [p[2] for p in solution.population], label="FSSH P2")
plot!(solution.t, [state for state in solution.state].-1, label="current_surface")

n_beads = 4
rpsh_sim = RingPolymerSimulation{FSSH}(atoms, Models.TullyModelOne(), n_beads; DoFs=1, temperature=10u"K")
r1 = RingPolymerArray(fill(-5.0, sim.DoFs, length(sim.atoms), n_beads))
v1 = RingPolymerArray(fill(8.9, sim.DoFs, length(sim.atoms), n_beads) ./ sim.atoms.masses[1])
rpsh_z = Dynamics.SurfaceHoppingVariables(v1, r1, 2, 1)

rpsh_solution = Dynamics.run_trajectory(rpsh_z, (0.0, 2500.0), rpsh_sim; output=(:population, :state))
plot!(rpsh_solution.t, [p[1] for p in rpsh_solution.population], label="RPSH P1")
plot!(rpsh_solution.t, [p[2] for p in rpsh_solution.population], label="RPSH P2")
plot!(rpsh_solution.t, [state for state in rpsh_solution.state].-1, label="current_surface")

nrpmd_sim = RingPolymerSimulation{NRPMD}(atoms, Models.TullyModelOne(), 1; DoFs=1, temperature=10u"K")
nrpmd_z = Dynamics.RingPolymerMappingDynamicals(v, r, 1, 2, 1)

nrpmd_solution = Dynamics.run_trajectory(nrpmd_z, (0.0, 2500.0), nrpmd_sim)

plot!(nrpmd_solution.t, hcat(Dynamics.get_population.(nrpmd_sim, nrpmd_solution.u)...)', label="nrpmd")

xlabel!("time")