using Test
using NonadiabaticMolecularDynamics

atoms = Atoms{Float64}([:C, :H])
cell = InfiniteCell{Float64}()

calc = Calculators.Free()

@test Simulation(3, 100, cell, atoms, calc, Dynamics.Classical()) isa Simulation
@test RingPolymerSimulation(3, 100, cell, atoms, calc, Dynamics.Classical(), 10) isa RingPolymerSimulation
