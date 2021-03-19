using Test
using NonadiabaticMolecularDynamics
using Unitful
using UnitfulAtomic

@test Dynamics.TwoTemperatureMDEF(x->exp(-x)) isa Dynamics.TwoTemperatureMDEF
atoms = Atoms([:H, :C])
sim = Simulation{MDEF}(atoms, Models.FrictionHarmonic(); temperature=10u"K", DoFs=1)

v = zeros(sim.DoFs, length(sim.atoms)) 
r = rand(sim.DoFs, length(sim.atoms)) 
u = ClassicalDynamicals(v, r)
du = zero(u)

sol = Dynamics.run_trajectory(u, (0.0, 100.0), sim; dt=1)
@test sol.u[1] ≈ u.x

f(t) = austrip(100u"K")*exp(-t)
two_temp = TwoTemperatureMDEF(f)
sim = Simulation(atoms, Models.FrictionHarmonic(), two_temp; DoFs=1)
sol = Dynamics.run_trajectory(u, (0.0, 100.0), sim; dt=1)
