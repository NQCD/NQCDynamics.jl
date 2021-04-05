using Test
using NonadiabaticMolecularDynamics

@testset "ClassicalDynamicals" begin
    v = rand(3, 10)
    r = rand(3, 10)
    z = ClassicalDynamicals(v, r)
    @test get_positions(z) == r # check the positions are extracted
    @test get_velocities(z) == v # check the momenta are extracted
    @test all(zero(z).x .== 0) # test zero function
    @test all(z.x .!= 0) # test that zero creates a copy and doesn't zero inplace
    
    @test get_flat_velocities(z) == get_velocities(z)[:]
    @test get_flat_positions(z) == get_positions(z)[:]
end
