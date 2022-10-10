using Test
using NQCDynamics
using Unitful
using UnitfulAtomic
using LinearAlgebra: diag
using ComponentArrays

atoms = Atoms([:H, :H])
model = CompositeFrictionModel(Free(), ConstantFriction(1, atoms.masses[1]))
sim = Simulation{MDEF}(atoms, model; temperature=10u"K")

v = zeros(size(sim))
r = randn(size(sim))
u = ComponentVector(v=v, r=r)
du = zero(u)

@testset "friction!" begin
    gtmp = zeros(length(r), length(r))
    NQCDynamics.DynamicsMethods.ClassicalMethods.friction!(gtmp, r, sim, 0.0)
    @test all(diag(gtmp) .≈ 1.0)
end

sol = run_dynamics(sim, (0.0, 100.0), u; output=OutputDynamicsVariables, dt=1)
@test sol[:OutputDynamicsVariables][1] ≈ u

f(t) = 100u"K"*exp(-ustrip(t))
model = CompositeFrictionModel(Harmonic(), RandomFriction(1))
sim = Simulation{MDEF}(atoms, model; temperature=f)
sol = run_dynamics(sim, (0.0, 100.0), u; output=OutputDynamicsVariables, dt=1)
 