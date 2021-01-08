using Test
using NonadiabaticMolecularDynamics

@testset "Phasespace" begin
    R = rand(3, 10)
    P = rand(3, 10)
    z = Phasespace(R, P)
    @test get_positions(z) == R # check the positions are extracted
    @test get_momenta(z) == P # check the momenta are extracted
    @test all(zero(z).x .== 0) # test zero function
    @test all(z.x .!= 0) # test that zero creates a copy and doesn't zero inplace
    
    @test get_flat_momenta(z) == get_momenta(z)[:]
    @test get_flat_positions(z) == get_positions(z)[:]
end

@testset "RingPolymerPhasespace" begin
    z = RingPolymerPhasespace(rand(3, 10), rand(3, 10), 20)
end