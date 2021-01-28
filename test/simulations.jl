using Test
using NonadiabaticMolecularDynamics

atoms = Atoms([:C, :H])
cell = InfiniteCell()

model = Models.Free()

@test Simulation(3, 100, cell, atoms, model, Classical()) isa Simulation
@test RingPolymerSimulation(3, 100, cell, atoms, model, Classical(), 10) isa RingPolymerSimulation

@test Simulation(atoms, model, Classical()) isa Simulation
@test RingPolymerSimulation(atoms, model, Classical(), 10) isa RingPolymerSimulation

@test Simulation(atoms, model) isa Simulation{Classical}
@test Simulation{Classical}(atoms, model) isa Simulation{Classical}
@test RingPolymerSimulation{Classical}(atoms, model, 10) isa RingPolymerSimulation{Classical}

# @test Simulation{MDEF}(atoms, model) isa Simulation{<:MDEF}

@test RingPolymerSimulation{NRPMD}(atoms, Models.DoubleWell(), 10) isa RingPolymerSimulation{<:NRPMD}
@test Simulation{Langevin}(atoms, model) isa Simulation{<:Langevin}
@test Simulation{FSSH}(atoms, Models.DoubleWell()) isa Simulation{<:FSSH}