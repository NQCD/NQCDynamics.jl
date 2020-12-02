using Test
using NonadiabaticMolecularDynamics
using DifferentialEquations
using Plots
using UnitfulAtomic
using PeriodicTable
using Random

atoms = Systems.AtomicParameters(Systems.InfiniteCell(), fill(:H, 10))
model = Models.Analytic.Harmonic(austrip(elements[:H].atomic_mass), 1e-3, 0)
p = Systems.System(atoms, model, 1)

z = Dynamics.Phasespace(randn(atoms.n_atoms), randn(atoms.n_atoms))

problem = ODEProblem(Dynamics.differential!, z, (0.0, 1e4), p)
solution = solve(problem)
plot(solution, vars=1:convert(Int, atoms.n_atoms))

