push!(LOAD_PATH, pwd())

using Test
using MDPrePostProcessing
using MDPrePostProcessing.Basic
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

@testset "System" begin
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
    @test get_positions(z) == R
    @test get_momenta(z) == P
end