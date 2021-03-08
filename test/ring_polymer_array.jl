
using NonadiabaticMolecularDynamics
using Test

b = RingPolymerArray([rand(3,2), rand(3,2), rand(3,2)])

@testset "transform!" begin
    b_original = deepcopy(b)
    transform!(b)
    @test !(b ≈ b_original)
    transform!(b)
    @test b ≈ b_original
end

@testset "get_centroid" begin
    b_original = deepcopy(b)
    centroid_original = get_centroid(b)
    transform!(b)
    centroid_middle = get_centroid(b)
    transform!(b)
    centroid_final = get_centroid(b)
    @test centroid_original ≈ centroid_middle ≈ centroid_final
end
