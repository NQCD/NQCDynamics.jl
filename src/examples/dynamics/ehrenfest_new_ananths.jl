using NQCDynamics
using Plots
#using Distributions

atoms = Atoms(2000)

# SINGLE TRAJECTORY
sim1 = Simulation{Ehrenfest}(atoms, AnanthModelOne())
sim2 = Simulation{Ehrenfest}(atoms, AnanthModelTwo())

k1 = 8.9
k2 = 16
k3 = 20
r = -5

u1 = DynamicsVariables(sim1, hcat(k1/atoms.masses[1]), hcat(r), SingleState(1, Adiabatic()))
u2 = DynamicsVariables(sim2, hcat(k2/atoms.masses[1]), hcat(r), SingleState(1, Adiabatic()))
output = Ensembles.OutputStateResolvedScattering1D(sim1, :adiabatic)
solution1 = run_trajectory(u1, (0.0, 2500.0), sim1; output=(output, :velocity))


new_population = [r.population[1] for r in solution1]
new_velocity = [r.velocity[1] for r in solution1]
print(solution1.transmission)
plot1 = plot(new_velocity, new_population)


solution2 = run_trajectory(u2, (0.0, 2500.0), sim2; output=(:population, :velocity))
plot2 = plot(solution2, :population)

plot(plot1, plot2, layout = (1, 2), size=(1000,400))
