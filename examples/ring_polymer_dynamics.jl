using Test
using NonadiabaticMolecularDynamics
using DifferentialEquations
using Plots
using UnitfulAtomic
using PeriodicTable
using Random
using Unitful

atoms = Systems.AtomicParameters(Systems.InfiniteCell(), [:H, :C])
model = Models.Analytic.Harmonic(austrip(elements[:H].atomic_mass), 1e-3, 10)
# model = Models.Analytic.Free()
beads = 10
DoF = 1
system = Systems.RingPolymerSystem(atoms, model, beads, DoF; temperature=100u"K")

R = cat([randn(DoF, 2) for i=1:beads]..., dims=3)
P = cat([zeros(DoF, 2) for i=1:beads]..., dims=3)
z = Dynamics.RingPolymerPhasespace(R, P)

problem = ODEProblem(Dynamics.differential!, z, (0.0, 1e4), system)
solution = solve(problem)
plot(solution)

