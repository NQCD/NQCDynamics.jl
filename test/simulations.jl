using Test
using NonadiabaticMolecularDynamics

atoms = Atoms([:C, :H])
cell = InfiniteCell()

model = Models.Free()

@test Simulation(3, 100, cell, atoms, model, Classical()) isa Simulation
@test RingPolymerSimulation(3, 100, cell, atoms, model, Classical(), 10) isa RingPolymerSimulation

@test Simulation(atoms, model, Classical()) isa Simulation
@test RingPolymerSimulation(atoms, model, Classical(), 10) isa RingPolymerSimulation
