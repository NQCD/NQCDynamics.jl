using Test
using NonadiabaticMolecularDynamics
using LinearAlgebra

@test_nowarn Systems.RingPolymerParameters{Float64}(10, 1, 1)
@test_nowarn Systems.RingPolymerParameters{Float64}(1, 1, 1)
@test_nowarn Systems.RingPolymerParameters{Float64}(2, 1, 1)
rp = Systems.RingPolymerParameters{Float64}(10, 1.0, vcat(fill(:C, 3), fill(:H, 5)), [:H])
@test rp.quantum_atoms == [4, 5, 6, 7, 8]
rp = Systems.RingPolymerParameters{Float64}(10, 1.0, 10)
@test rp.quantum_atoms == collect(1:10)

@test sort(eigvals(rp.springs)) ≈ sort(rp.normal_mode_springs)

atoms = Systems.AtomicParameters(Systems.InfiniteCell(), vcat(fill(:H, 10), :O))
Systems.RingPolymerSystem(atoms, Models.Analytic.Free(), 10, 10)
system = Systems.RingPolymerSystem(atoms, Models.Analytic.Free(), 10, 10)
@test system.ring_polymer.quantum_atoms == collect(1:11)
system = Systems.RingPolymerSystem(atoms, Models.Analytic.Free(), 10, 10, quantum_nuclei=[:H])
@test system.ring_polymer.quantum_atoms == collect(1:10)

