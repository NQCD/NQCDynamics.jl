using Test
using NonadiabaticMolecularDynamics
using DifferentialEquations
using Plots
using UnitfulAtomic
using PeriodicTable
using Random

atoms = Systems.AtomicParameters(Systems.InfiniteCell(), fill(:H, 10))
model = Models.Analytic.Harmonic(austrip(elements[:H].atomic_mass), 1e-1, 0)
p = Systems.System{Langevin}(atoms, model, 1u"K", 1, 1)

z = Dynamics.Phasespace(randn(n_atoms(p)), randn(n_atoms(p)))

problem = SDEProblem(Dynamics.differential!, Dynamics.random_force!, z, (0.0, 1e4), p)
solution = solve(problem)
plot(solution)

