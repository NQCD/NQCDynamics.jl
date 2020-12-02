using Test
using NonadiabaticMolecularDynamics
using DifferentialEquations
using Plots
using UnitfulAtomic
using PeriodicTable
using Random

atoms = Systems.AtomicParameters(Systems.InfiniteCell(), fill(:H, 2))
model = Models.Analytic.Harmonic(austrip(elements[:H].atomic_mass), 1e-3, 10)
# model = Models.Analytic.Free()
n_beads = 10
temperature = 1
system = Systems.RingPolymerSystem(atoms, model, 9, 1e-3, 1)

R = randn(system.atomic_parameters.n_atoms * system.ring_polymer.n_beads)
P = randn(system.atomic_parameters.n_atoms * system.ring_polymer.n_beads).*10
z = Dynamics.Phasespace(R, P)

problem = ODEProblem(Dynamics.differential!, z, (0.0, 1e4), system)
solution = solve(problem)
plot(solution, vars=1:convert(Int, system.atomic_parameters.n_atoms*system.ring_polymer.n_beads))

