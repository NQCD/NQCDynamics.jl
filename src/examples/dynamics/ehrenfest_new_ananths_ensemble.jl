using NQCDynamics


atoms = Atoms(1980)
sim1 = Simulation{Ehrenfest}(atoms, AnanthModelOne())

e = 0.03
k = sqrt(e*2*atoms.masses[1])
r = Normal(-5, 1/sqrt(0.25))
v = k / atoms.masses[1]
distribution = DynamicalDistribution(v, r, size(sim1))* SingleState(1, Adiabatic())


n_traj = 500
tspan = (0.0, 3000.0)
solution1 = run_ensemble(sim1, tspan, distribution; 
    trajectories=n_traj, output=:velocity)
new_velocity = [r.velocity[end] for r in solution1]
results = reduce(vcat, new_velocity*atoms.masses[1])

#using Plots
using StatsPlots
plot(density(results))
xlims!(-20,20)