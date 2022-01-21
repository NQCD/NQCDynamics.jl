using NonadiabaticMolecularDynamics
using Test
using StatsBase: mean

b = RingPolymerArray(rand(2,3,4), quantum=[2])
beads = RingPolymers.RingPolymerParameters{Float64}(4, 1.0, [2])

@testset "transform!" begin
    b_original = deepcopy(b)
    RingPolymers.transform!(beads, b)
    @test !(b ≈ b_original)
    RingPolymers.transform!(beads, b)
    @test b ≈ b_original
end

@testset "setindex!" begin
    b[1,1,1] = 1.0
    @test all(b[1,1,:] .≈ 1.0) # Test all classical atoms set to 1
    b[1,2,1] = 1.0
    !(b[1,2,2] ≈ 1.0) # Test quantum atoms not all 1
end

@testset "similar" begin
    c = similar(b)
    @test c.quantum_atoms == b.quantum_atoms
    @test c.classical_atoms == b.classical_atoms
    @test c.normal == b.normal
end

@testset "get_centroid" begin
    b_original = deepcopy(b)
    centroid_original = RingPolymers.get_centroid(b)
    RingPolymers.transform!(beads, b)
    @test !(b ≈ b_original)
    centroid_middle = RingPolymers.get_centroid(b)
    RingPolymers.transform!(beads, b)
    @test b ≈ b_original
    centroid_final = RingPolymers.get_centroid(b)

    @test centroid_original ≈ centroid_middle ≈ centroid_final

    R = rand(size(b)...)
    centroid = dropdims(mean(R; dims=3); dims=3)
    @test centroid ≈ RingPolymers.get_centroid(R)
    out = similar(centroid) .+ 10
    RingPolymers.get_centroid!(out, R)
    @test centroid ≈ out
end

@testset "broadcasting" begin
    data = rand(2,3,4)
    A = RingPolymerArray(data, quantum=[2])

    test_data = copy(data)
    RingPolymers.constrain_classical_atoms!(test_data, A.classical_atoms)

    B = A ./ 10
    @test test_data ./ 10 ≈ B.data

    A ./= 10
    @test test_data ./ 10 ≈ A.data
end
