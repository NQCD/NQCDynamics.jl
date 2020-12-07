using Test
using NonadiabaticMolecularDynamics
using DifferentialEquations
using Plots
using UnitfulAtomic
using PeriodicTable
using Random
using Unitful

atoms = Atoms.AtomicParameters(Atoms.InfiniteCell(), fill(:H, 1))
model = Models.Analytic.Harmonic(austrip(elements[:H].atomic_mass), 1, 0)
p = System{Langevin}(atoms, model, 300u"K", 0.1, 1)

z = Dynamics.Phasespace(randn(n_atoms(p)), randn(n_atoms(p)))

problem = SDEProblem(Dynamics.differential!, Dynamics.random_force!, z, (0.0, 50), p)
solution = solve(problem, EM(), dt=0.001)
plot(solution)

