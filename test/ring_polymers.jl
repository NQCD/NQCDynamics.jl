using Test
using NonadiabaticMolecularDynamics
using LinearAlgebra

@test_nowarn RingPolymers.RingPolymerParameters{Float64}(10, 1, 1)
@test_nowarn RingPolymers.RingPolymerParameters{Float64}(1, 1, 1)
@test_nowarn RingPolymers.RingPolymerParameters{Float64}(2, 1, 1)
rp = RingPolymers.RingPolymerParameters{Float64}(10, 1.0, vcat(fill(:C, 3), fill(:H, 5)), [:H])
@test rp.quantum_atoms == [4, 5, 6, 7, 8]
rp = RingPolymers.RingPolymerParameters{Float64}(10, 1.0, 10)
@test rp.quantum_atoms == collect(1:10)

@test sort(eigvals(rp.springs)) ≈ sort(rp.normal_mode_springs)

atoms = Atoms{Float64}(vcat(fill(:H, 10), :O))
model = NonadiabaticModels.Free(3)
sim = RingPolymerSimulation{Classical}(atoms, model, 10)
@test sim.beads.quantum_atoms == collect(1:11)
sim = RingPolymerSimulation{Classical}(atoms, model, 10; quantum_nuclei=[:H])
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
RingPolymers.transform_to_normal_modes!(sim.beads, R)
@test !(R ≈ R_original)
RingPolymers.transform_from_normal_modes!(sim.beads, R)
@test R ≈ R_original

half = RingPolymers.cayley_propagator(sim.beads, 0.1; half=true)
full = RingPolymers.cayley_propagator(sim.beads, 0.1; half=false)
@test full ≈ half .* half
