using Test
using NonadiabaticMolecularDynamics

atoms = Atoms{Float64}([:C, :H])
cell = InfiniteCell{Float64}()

model = Models.Free()

@test Simulation(3, 100, cell, atoms, model, Dynamics.Classical()) isa Simulation
@test RingPolymerSimulation(3, 100, cell, atoms, model, Dynamics.Classical(), 10) isa RingPolymerSimulation
