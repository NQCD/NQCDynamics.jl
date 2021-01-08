using NonadiabaticMolecularDynamics
using Plots

atoms = Atoms{Float64}([:H])
sim = Simulation(atoms, Models.TullyModelOne(), Dynamics.FSSH{Float64}(1, 1, 2); DoFs=1)

R = fill(-5.0, sim.DoFs, length(sim.atoms)) 
P = fill(8.9, sim.DoFs, length(sim.atoms)) 
z = SurfaceHoppingPhasespace(R, P, 2, 1)

solution = Dynamics.run_trajectory(z, (0.0, 2500.0), sim)

plot(solution.t, [real(Dynamics.get_density_matrix(u)[1,1]) for u in solution.u], label="σ[1,1]")
plot!(solution.t, [real(Dynamics.get_density_matrix(u)[2,2]) for u in solution.u], label="σ[2,2]")
plot!(solution.t, [u.state for u in solution.u].-1, label="current surface")

xlabel!("time")