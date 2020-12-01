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
    p = AtomicParameters(cell, [:O, :C, :H, :Pd])
    @test p.n_atoms == 4
    @test p.atom_types[1] == :O
    @test elements[p.atom_types[1]] == elements[:O]
end

@testset "System" begin
    p = AtomicParameters(Cell(zeros(3, 3)), [:C, :C])
    system = System(p, Free())
end