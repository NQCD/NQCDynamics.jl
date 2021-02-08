using Test
using NonadiabaticMolecularDynamics
using Unitful

@test Dynamics.MDEF{Float64}(10) isa Dynamics.MDEF
@test Dynamics.TwoTemperatureMDEF{Float64}(10, x->exp(-x)) isa Dynamics.TwoTemperatureMDEF
atoms = Atoms([:H, :C])
sim = Simulation{MDEF}(atoms, Models.FrictionHarmonic(); temperature=10u"K", DoFs=1)

v = zeros(sim.DoFs, length(sim.atoms)) 
r = rand(sim.DoFs, length(sim.atoms)) 
u = ClassicalDynamicals(v, r)
du = zero(u)

sol = Dynamics.run_trajectory(u, (0.0, 100.0), sim; dt=1)
@test sol.u[1] ≈ u.x
