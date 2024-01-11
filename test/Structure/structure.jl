using Test
using NQCDynamics
using Unitful,UnitfulAtomic
using LinearAlgebra

@testset "Minimum distance translation" begin
    cell=PeriodicCell(10 .* I(3)) # 10x10x10 cell
    positions=[0 0; 0 0; 0 49]

    positions_new=copy(positions)
    positions_new[:,2]+=Structure.minimum_distance_translation(positions, 1, 2, cell)

    @test positions_new[:,2] == [0.0,0.0,-1.0]
end

@testset "PBC-compatible functions" begin
    cell=PeriodicCell(10 .* I(3)) # 10x10x10 cell
    atoms=Atoms([:H, :H])

    positions=[0 0; 0 0; 0 9]
    @test 1.0 == austrip(Structure.pbc_distance(positions, 1, 2, cell))
    sim=Simulation(atoms, Free(3); cell=cell)
    @test 1.0 == austrip(Structure.pbc_distance(positions, 1,2,sim))
    positions=[0 0; 0 0; 1 9]
    @test Structure.pbc_center_of_mass(positions, 1, 2, cell) == ([0.0, 0.0, 0.0] || [0.0 ,0.0 ,10.0])
    @test Structure.pbc_center_of_mass(positions, 1, 2, sim) == ([0.0, 0.0, 0.0] || [0.0 ,0.0 ,10.0])
end