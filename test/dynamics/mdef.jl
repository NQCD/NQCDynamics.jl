using Test
using NonadiabaticMolecularDynamics
using Unitful
using UnitfulAtomic
using RecursiveArrayTools
using LinearAlgebra: diag

atoms = Atoms([:H, :H])
@test Dynamics.TwoTemperatureMDEF(atoms.masses, 3, x->exp(-x)) isa Dynamics.TwoTemperatureMDEF
sim = Simulation{MDEF}(atoms, Models.FreeConstantFriction(atoms.masses[1]); temperature=10u"K", DoFs=2)

v = zeros(sim.DoFs, length(sim.atoms))
r = rand(sim.DoFs, length(sim.atoms))
u = ClassicalDynamicals(v, r)
du = zero(u)

@testset "friction!" begin
    gtmp = ArrayPartition(zeros(length(r), length(r)), zeros(length(r)))
    NonadiabaticMolecularDynamics.Dynamics.friction!(gtmp, r, sim, 0.0)
    @test all(diag(gtmp.x[1]) .≈ 1.0)
end

sol = Dynamics.run_trajectory(u, (0.0, 100.0), sim; dt=1)
@test sol.u[1] ≈ u.x

f(t) = austrip(100u"K")*exp(-t)
sim = Simulation{TwoTemperatureMDEF}(atoms, Models.FrictionHarmonic(), f; DoFs=2)
sol = Dynamics.run_trajectory(u, (0.0, 100.0), sim; dt=1)
