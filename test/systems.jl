using Test
using NonadiabaticMolecularDynamics.Systems
using NonadiabaticMolecularDynamics.Atoms
using NonadiabaticMolecularDynamics.Models
using Unitful
using UnitfulAtomic
using PeriodicTable

@testset "Cell" begin
    x = [1.0, 0.0, 0.0]
    y = [0.0, 1.0, 0.0]
    z = [0.0, 0.0, 1.0]
    @test PeriodicCell([x y z] .* u"Å", [true, true, true]) isa AbstractCell
    @test PeriodicCell([x y z], u"bohr").vectors == PeriodicCell([x y z]).vectors # Defaults to Angstrom
    @test PeriodicCell([x y z], u"m") isa AbstractCell
    @test_throws MethodError PeriodicCell([x y z], u"s") # Throws error for non-length units
    
    a = PeriodicCell([x y z])
    @test a.periodicity == [true, true, true]
    set_periodicity!(a, [false, true, false])
    @test a.periodicity == [false, true, false]
end

@testset "Parameters" begin
    x = [1, 0, 0]
    y = [0, 1, 0]
    z = [0, 0, 1]
    cell = PeriodicCell([x y z])
    p = AtomicParameters(cell, [:O, :C, :H, :Pd])
    @test p.n_atoms == 4
    @test p.atom_types[1] == :O
    @test elements[p.atom_types[1]] == elements[:O]
end

@testset "System" begin
    p = AtomicParameters(PeriodicCell(zeros(3, 3)), [:C, :C])
    system = System(p, Models.Analytic.Free())
end