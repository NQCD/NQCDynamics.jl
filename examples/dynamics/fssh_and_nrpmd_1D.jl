using NonadiabaticMolecularDynamics
using Plots

atoms = Atoms([:H])
sim = Simulation{FSSH}(atoms, Models.TullyModelOne(); DoFs=1)

r = fill(-5.0, sim.DoFs, length(sim.atoms))
v = fill(8.9, sim.DoFs, length(sim.atoms)) ./ sim.atoms.masses[1]
z = SurfaceHoppingDynamicals(v, r, 2, 1)

solution = Dynamics.run_trajectory(z, (0.0, 2500.0), sim)

plot(solution.t, [real(get_density_matrix(u)[1,1]) for u in solution.u], label="σ[1,1]")
plot!(solution.t, [real(get_density_matrix(u)[2,2]) for u in solution.u], label="σ[2,2]")
plot!(solution.t, [u.state for u in solution.u].-1, label="current surface")

nrpmd_sim = RingPolymerSimulation{NRPMD}(atoms, Models.TullyModelOne(), 1; DoFs=1)
nrpmd_z = Dynamics.RingPolymerMappingDynamicals(v, r, 1, 2, 1)

nrpmd_solution = Dynamics.run_trajectory(nrpmd_z, (0.0, 2500.0), nrpmd_sim)

plot!(nrpmd_solution.t, hcat(Dynamics.get_population.(nrpmd_solution.u)...)', label="nrpmd")

xlabel!("time")