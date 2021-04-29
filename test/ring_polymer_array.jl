using NonadiabaticMolecularDynamics
using Test

b = RingPolymerArray(rand(2,3,4), quantum=[2])
beads = RingPolymerParameters{Float64}(4, 1.0, [2])

@testset "transform!" begin
    b_original = deepcopy(b)
    transform!(beads, b)
    @test !(b ≈ b_original)
    transform!(beads, b)
    @test b ≈ b_original
end

@testset "setindex!" begin
    b[1,1,1] = 1.0
    @test all(b[1,1,:] .≈ 1.0) # Test all classical atoms set to 1
    b[1,2,1] = 1.0
    !(b[1,2,2] ≈ 1.0) # Test quantum atoms not all 1
end

@testset "get_centroid" begin
    b_original = deepcopy(b)
    centroid_original = get_centroid(b)
    transform!(beads, b)
    @test !(b ≈ b_original)
    centroid_middle = get_centroid(b)
    transform!(beads, b)
    @test b ≈ b_original
    centroid_final = get_centroid(b)

    @test centroid_original ≈ centroid_middle ≈ centroid_final
end
