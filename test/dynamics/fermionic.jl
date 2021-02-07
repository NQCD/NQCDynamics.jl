using Test
using NonadiabaticMolecularDynamics
using Unitful

f = Dynamics.FermionicBath(1.0, 1.0, 1.0)

sim = RingPolymerSimulation(Atoms{Float64}([:H]), Models.Free(), f, 10; temperature=10u"K")

R = rand(3, 1, 10)
evaluate_potential_energy(sim, R)