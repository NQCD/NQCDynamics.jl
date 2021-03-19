using Test
using NonadiabaticMolecularDynamics
using LinearAlgebra

@test_nowarn RingPolymerParameters{Float64}(10, 1, 1)
@test_nowarn RingPolymerParameters{Float64}(1, 1, 1)
@test_nowarn RingPolymerParameters{Float64}(2, 1, 1)
rp = RingPolymerParameters{Float64}(10, 1.0, vcat(fill(:C, 3), fill(:H, 5)), [:H])
@test rp.quantum_atoms == [4, 5, 6, 7, 8]
rp = RingPolymerParameters{Float64}(10, 1.0, 10)
@test rp.quantum_atoms == collect(1:10)

@test sort(eigvals(rp.springs)) ≈ sort(rp.normal_mode_springs)

atoms = Atoms{Float64}(vcat(fill(:H, 10), :O))
cell = InfiniteCell()
model = Models.Free()
sim = RingPolymerSimulation(3, 1.0, cell, atoms, model, Dynamics.Classical(), 10)
@test sim.beads.quantum_atoms == collect(1:11)
sim = RingPolymerSimulation(3, 1.0, cell, atoms, model, Dynamics.Classical(), 10, [:H])
@test sim.beads.quantum_atoms == collect(1:10)

U = sim.beads.U
S = sim.beads.springs
λ = sim.beads.normal_mode_springs
@test U'U ≈ I
@test abs(det(U)) ≈ 1
@test U'S*U ≈ diagm(λ)

R = rand(3, 11, 10)
R_original = copy(R)
@test R == R_original
transform_to_normal_modes!(sim.beads, R)
@test !(R ≈ R_original)
transform_from_normal_modes!(sim.beads, R)
@test R ≈ R_original