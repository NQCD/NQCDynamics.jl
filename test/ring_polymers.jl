using Test
using NonadiabaticMolecularDynamics
using LinearAlgebra

@test_nowarn Systems.RingPolymerParameters{Float64}(10, 1)
@test_nowarn Systems.RingPolymerParameters{Float64}(1, 1)
@test_nowarn Systems.RingPolymerParameters{Float64}(2, 1)

rp = Systems.RingPolymerParameters(10, 1.0)
@test sort(eigvals(rp.springs)) ≈ sort(rp.normal_mode_springs)

atoms = Systems.AtomicParameters(Systems.InfiniteCell(), fill(:H, 10))
@test_nowarn Systems.RingPolymerSystem(atoms, Free(), 10, 10)
system = Systems.RingPolymerSystem(atoms, Free(), 10, 10)