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
system = Systems.RingPolymerSystem(atoms, model, n_beads, 5e-4, 1)

R = randn(Systems.n_atoms(system), n_beads)
P = zeros(Systems.n_atoms(system), n_beads)
z = Dynamics.RingPolymerPhasespace(R, P)

problem = ODEProblem(Dynamics.differential!, z, (0.0, 1e4), system)
solution = solve(problem)
plot(solution, vars=1:convert(Int, 2Systems.n_atoms(system)*n_beads))

