using Test
using NonadiabaticMolecularDynamics

@test Dynamics.NRPMD{Float64}(10) isa Dynamics.NRPMD
atoms = Atoms([:H])
sim = RingPolymerSimulation{NRPMD}(atoms, Models.DoubleWell(), 10; DoFs=1)

v = RingPolymerArray(zeros(sim.DoFs, length(sim.atoms), length(sim.beads)))
r = RingPolymerArray(rand(sim.DoFs, length(sim.atoms), length(sim.beads)))
u = Dynamics.RingPolymerMappingVariables(v, r, 2, 2)
qmap = Dynamics.get_mapping_positions(u)
pmap = Dynamics.get_mapping_momenta(u)
population = Dynamics.get_population(sim, u)
@test population[1] ≈ 0 atol=1e-10
@test population[2] ≈ 1

du = zero(u)
Dynamics.motion!(du, u, sim, 0.0)

Dynamics.run_trajectory(u, (0, 1.0), sim)
