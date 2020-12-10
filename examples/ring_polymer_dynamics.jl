using Test
using NonadiabaticMolecularDynamics
using DifferentialEquations
using Plots
using UnitfulAtomic
using PeriodicTable
using Random

atoms = Systems.AtomicParameters(Systems.InfiniteCell(), [:H, :C])
model = Models.Analytic.Harmonic(austrip(elements[:H].atomic_mass), 1e-3, 10)
# model = Models.Analytic.Free()
n_beads = 10
temperature = 1
n_DoF = 1
system = Systems.RingPolymerSystem(atoms, model, n_beads, 5e-4, n_DoF)

R = ArrayPartition([randn(n_DoF, 2) for i=1:n_beads]...)
P = ArrayPartition([zeros(n_DoF, 2) for i=1:n_beads]...)
z = Dynamics.RingPolymerPhasespace(R, P)

problem = ODEProblem(Dynamics.differential!, z, (0.0, 1e4), system)
solution = solve(problem)
plot(solution)

