using Test
using NonadiabaticMolecularDynamics.Systems
using Unitful
using UnitfulAtomic
using PeriodicTable

@testset "Cell" begin
    x = [1.0, 0.0, 0.0]
    y = [0.0, 1.0, 0.0]
    z = [0.0, 0.0, 1.0]
    @test Cell([x y z] .* u"Å") isa Cell
    @test Cell([x y z], u"bohr").vectors == Cell([x y z]).vectors # Defaults to Angstrom
    @test Cell([x y z], u"m") isa Cell
    @test_throws MethodError Cell([x y z], u"s") # Throws error for non-length units
end

@testset "Parameters" begin
    x = [1, 0, 0]
    y = [0, 1, 0]
    z = [0, 0, 1]
    cell = Cell([x y z])
    p = SystemParameters(cell, [:O, :C, :H, :Pd])
    @test p.n_beads == 1
    @test p.n_atoms == 4
    @test p.atom_types[1] == :O
    @test elements[p.atom_types[1]] == elements[:O]
end

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

@testset "System" begin
    p = SystemParameters(Cell(zeros(3, 3)), [:C, :C])
    z = Phasespace(rand(p.n_atoms), rand(p.n_atoms))
    system = System(p, z)
    @test system.parameters == p
    @test system.dynamical_variables == z
end