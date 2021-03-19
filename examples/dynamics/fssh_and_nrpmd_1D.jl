using NonadiabaticMolecularDynamics
using Plots
using Unitful

atoms = Atoms([:H])
sim = Simulation{FSSH}(atoms, Models.TullyModelOne(); DoFs=1)

r = fill(-5.0, sim.DoFs, length(sim.atoms))
v = fill(8.9, sim.DoFs, length(sim.atoms)) ./ sim.atoms.masses[1]
z = SurfaceHoppingDynamicals(v, r, 2, 1)

solution = Dynamics.run_trajectory(z, (0.0, 2500.0), sim; output=(:density_matrix, :state))

plot(solution.t, [real(σ[1,1]) for σ in solution.density_matrix], label="σ[1,1]")
plot!(solution.t, [real(σ[2,2]) for σ in solution.density_matrix], label="σ[2,2]")
plot!(solution.t, [state for state in solution.state].-1, label="current_surface")

n_beads = 4
rpsh_sim = RingPolymerSimulation{FSSH}(atoms, Models.TullyModelOne(), n_beads; DoFs=1, temperature=10u"K")
r1 = fill(-5.0, sim.DoFs, length(sim.atoms), n_beads)
v1 = fill(8.9, sim.DoFs, length(sim.atoms), n_beads) ./ sim.atoms.masses[1]
rpsh_z = Dynamics.SurfaceHoppingDynamicals(v1, r1, 2, 1)

rpsh_solution = Dynamics.run_trajectory(rpsh_z, (0.0, 2500.0), rpsh_sim; output=(:density_matrix, :state))

plot(rpsh_solution.t, [real(σ[1,1]) for σ in rpsh_solution.density_matrix], label="σ[1,1]")
plot!(rpsh_solution.t, [real(σ[2,2]) for σ in rpsh_solution.density_matrix], label="σ[2,2]")
plot!(rpsh_solution.t, [state for state in rpsh_solution.state].-1, label="current_surface")

nrpmd_sim = RingPolymerSimulation{NRPMD}(atoms, Models.TullyModelOne(), 1; DoFs=1, temperature=10u"K")
nrpmd_z = Dynamics.RingPolymerMappingDynamicals(v, r, 1, 2, 1)

nrpmd_solution = Dynamics.run_trajectory(nrpmd_z, (0.0, 2500.0), nrpmd_sim)

plot!(nrpmd_solution.t, hcat(Dynamics.get_population.(nrpmd_solution.u)...)', label="nrpmd")

xlabel!("time")