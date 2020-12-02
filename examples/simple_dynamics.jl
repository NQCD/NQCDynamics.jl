using Test
using NonadiabaticMolecularDynamics
using DifferentialEquations
using Plots
using UnitfulAtomic
using PeriodicTable
using Random

x = [1, 0, 0]
y = [0, 1, 0]
z = [0, 0, 1]
atoms = Systems.AtomicParameters(Systems.PeriodicCell([x y z]), fill(:H, 10))
model = Models.Analytic.Harmonic(austrip(elements[:H].atomic_mass), 1e-3, 0)
p = Systems.System(atoms, model, 1)

z = Dynamics.Phasespace(randn(atoms.n_atoms), randn(atoms.n_atoms))
du = zero(z)
Dynamics.differential!(du, z, p, 0)

problem = ODEProblem(Dynamics.differential!, z, (0.0, 1e4), p)
solution = solve(problem)
plot(solution, vars=1:convert(Int, atoms.n_atoms))

