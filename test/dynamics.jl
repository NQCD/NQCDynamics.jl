using Test
using NonadiabaticMolecularDynamics.Dynamics
using Unitful

@testset "Phasespace" begin
    R = rand(10)
    P = rand(10)
    z = Phasespace(R, P, (u"m", u"kg*m/s"))
    @test get_positions(z) != R
    @test get_momenta(z) != R
    z = Phasespace(R, P)
    @test get_positions(z) == R # check the positions are extracted
    @test get_momenta(z) == P # check the momenta are extracted
    @test all(zero(z).x .== 0) # test zero function
    @test all(z.x .!= 0) # test that zero creates a copy and doesn't zero inplace
end

@testset "MappingPhasespace" begin
    R = rand(10)
    P = rand(10)
    z = MappingPhasespace(R, P, 2, 1)
    
    @test get_positions(z) == R
    @test get_momenta(z) == P
    get_mapping_positions(z)
    get_mapping_momenta(z)
end