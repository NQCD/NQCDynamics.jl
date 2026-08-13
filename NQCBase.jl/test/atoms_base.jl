using NQCBase
using AtomsBase: AtomsBase
using AtomsBaseTesting: AtomsBaseTesting
using Unitful, UnitfulAtomic
using Test

hydrogen = AtomsBase.isolated_system([:H => [0, 0, 1.]u"bohr",
                            :H => [0, 0, 3.]u"bohr"])

box = 10.26 / 2 * [[0, 0, 1], [1, 0, 1], [1, 1, 0]]u"bohr"
silicon = AtomsBase.periodic_system([:Si =>  ones(3)/8,
                           :Si => -ones(3)/8],
                           box, fractional=true)

@testset "Atoms conversions" begin
    @test Atoms(hydrogen) == Atoms([:H, :H])
    @test Atoms(silicon) == Atoms([:Si, :Si])
end

@testset "Cell conversions" begin
    @test Cell(hydrogen) === InfiniteCell()
    @test Cell(silicon) isa PeriodicCell
end

@testset "System" begin
    @test System(Atoms([:H, :H]), rand(3,2)) isa AtomsBase.FlexibleSystem
    @test System(Atoms([:H, :H]), rand(3,2), rand(3,2)) isa AtomsBase.FlexibleSystem
    @test System(Atoms([:H, :H]), rand(3,2), rand(3,2), Cell(silicon)) isa AtomsBase.FlexibleSystem
end

@testset "Forward and backward conversion" begin
    atoms = Atoms(silicon)
    cell = Cell(silicon)
    position = Position(silicon)
    velocity = Velocity(silicon)
    AtomsBaseTesting.test_approx_eq(System(atoms, position, cell), silicon)
    AtomsBaseTesting.test_approx_eq(System(atoms, position, velocity, cell), silicon)
end

@testset "Trajectory" begin
    position = [rand(3,3) for i in 1:10]
    velocity = [rand(3,3) for i in 1:10]
    Trajectory(Atoms([:H, :C, :N]), position, velocity) isa Vector{<:AtomsBase.FlexibleSystem}
    Trajectory(Atoms([:H, :C, :N]), position, velocity, Cell(silicon)) isa Vector{<:AtomsBase.FlexibleSystem}
end
