using Test
using NonadiabaticMolecularDynamics.Systems
using NonadiabaticMolecularDynamics.Dynamics
using NonadiabaticMolecularDynamics.Electronics
using NonadiabaticMolecularDynamics.Models
using DifferentialEquations
using Plots
using UnitfulAtomic
using PeriodicTable

x = [1, 0, 0]
y = [0, 1, 0]
z = [0, 0, 1]
system = SystemParameters(Cell([x y z]), fill(:H, 10))
z = Phasespace(randn(system.n_atoms), randn(system.n_atoms))

model = Models.Analytic.Harmonic(austrip(elements[:H].atomic_mass), 1, 0)
p = Dynamics.DynamicsParameters(system, model)

du = zero(z)
Dynamics.differential!(du, z, p, 0)

problem = ODEProblem(Dynamics.differential!, z, (0.0, 10.0), p)
solution = solve(problem)
@info get_positions(z)
plot(solution, vars=1:convert(Int, system.n_atoms))
