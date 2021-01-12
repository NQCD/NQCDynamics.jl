using Test
using NonadiabaticMolecularDynamics

@test Dynamics.NRPMD{Float64}(10) isa Dynamics.NRPMD
atoms = Atoms{Float64}([:H])
sim = RingPolymerSimulation(atoms, Models.DoubleWell(), Dynamics.NRPMD{Float64}(2), 10; DoFs=1)

R = zeros(sim.DoFs, length(sim.atoms), length(sim.beads)) 
P = rand(sim.DoFs, length(sim.atoms), length(sim.beads)) 
u = Dynamics.RingPolymerMappingPhasespace(R, P, 2, 2)
qmap = Dynamics.get_mapping_positions(u)
pmap = Dynamics.get_mapping_momenta(u)
population = Dynamics.get_population(u)
@test population[1] ≈ 0 atol=1e-10
@test population[2] ≈ 1

du = zero(u)
Dynamics.motion!(du, u, sim, 0.0)

Dynamics.run_trajectory(u, (0, 1.0), sim)
