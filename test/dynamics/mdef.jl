using Test
using NonadiabaticMolecularDynamics
using Unitful
using UnitfulAtomic
using RecursiveArrayTools: ArrayPartition
using LinearAlgebra: diag

atoms = Atoms([:H, :H])
sim = Simulation{MDEF}(atoms, ConstantFriction(Free(),atoms.masses[1]); temperature=10u"K", DoFs=2)

v = zeros(sim.DoFs, length(sim.atoms))
r = rand(sim.DoFs, length(sim.atoms))
u = ComponentVector(v=v, r=r)
du = zero(u)

@testset "friction!" begin
    gtmp = zeros(length(r), length(r))
    NonadiabaticMolecularDynamics.Dynamics.ClassicalMethods.friction!(gtmp, r, sim, 0.0)
    @test all(diag(gtmp) .≈ 1.0)
end

sol = Dynamics.run_trajectory(u, (0.0, 100.0), sim; dt=1)
@test sol.u[1] ≈ u

f(t) = 100u"K"*exp(-ustrip(t))
sim = Simulation{MDEF}(atoms, NonadiabaticModels.RandomFriction(Harmonic()); DoFs=2, temperature=f)
sol = Dynamics.run_trajectory(u, (0.0, 100.0), sim; dt=1)
