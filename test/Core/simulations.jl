using Test
using NQCDynamics
using Unitful
using Distributions
using PyCall

atoms = Atoms([:C, :H])
cell = InfiniteCell()

model = NQCModels.Free()

@test RingPolymerSimulation(100, cell, atoms, model, Classical(), 10, [:H]) isa RingPolymerSimulation

@test Simulation(atoms, model, Classical()) isa Simulation
@test RingPolymerSimulation(atoms, model, Classical(), 10) isa RingPolymerSimulation

@test Simulation(atoms, model) isa Simulation{Classical}
@test Simulation{Classical}(atoms, model) isa Simulation{Classical}
@test RingPolymerSimulation{Classical}(atoms, model, 10) isa RingPolymerSimulation{Classical}

# @test Simulation{MDEF}(atoms, model) isa Simulation{<:MDEF}

@test RingPolymerSimulation{NRPMD}(atoms, NQCModels.DoubleWell(), 10) isa RingPolymerSimulation{<:NRPMD}
@test Simulation{Langevin}(atoms, model) isa Simulation{<:Langevin}
@test Simulation{FSSH}(atoms, NQCModels.DoubleWell()) isa Simulation{<:FSSH}

@testset "get_temperature" begin
    sim = Simulation(atoms, model; temperature=10u"K")
    @test NQCDynamics.get_temperature(sim) isa Real
    sim = Simulation(atoms, model; temperature=x->10u"K")
    @test NQCDynamics.get_temperature(sim) isa Real
    sim = Simulation(atoms, model; temperature=10)
    @test NQCDynamics.get_temperature(sim, 1) isa Real
end

@testset "get_ring_polymer_temperature" begin
    sim = RingPolymerSimulation(atoms, model, 10; temperature=10u"K")
    @test NQCDynamics.get_ring_polymer_temperature(sim) isa Real
    sim = RingPolymerSimulation(atoms, model, 10; temperature=x->10u"K")
    @test NQCDynamics.get_ring_polymer_temperature(sim) isa Real
    sim = RingPolymerSimulation(atoms, model, 10; temperature=10)
    @test NQCDynamics.get_ring_polymer_temperature(sim, 1) isa Real
end

@testset "nfunctions and size" begin
    sim = Simulation(atoms, model)
    @test natoms(sim) == 2
    @test ndofs(sim) == 1
    @test_throws MethodError nbeads(sim)
    @test size(sim) == (1, 2)

    sim = RingPolymerSimulation(atoms, model, 10)
    @test natoms(sim) == 2
    @test ndofs(sim) == 1
    @test nbeads(sim) == 10
    @test size(sim) == (1, 2, 10)
end

@testset "masses" begin
    sim = Simulation(atoms, model)
    @test masses(sim) == atoms.masses
    @test masses(sim, 2) == atoms.masses[2]
    @test masses(sim, CartesianIndex(1, 2)) == atoms.masses[2]

    sim = RingPolymerSimulation(atoms, model, 5)
    @test masses(sim) == atoms.masses
    @test masses(sim, 2) == atoms.masses[2]
    @test masses(sim, CartesianIndex(1, 2, 3)) == atoms.masses[2]
end

@testset "mobileatoms and distribution generation" begin
    sim = Simulation(atoms, model)
    @test NQCModels.mobileatoms(sim) == 1:2
    @test NQCDistributions.VelocityBoltzmann(10, sim).sampleable == VelocityBoltzmann(10, atoms.masses, size(sim)).sampleable

    # Generate a structure with constraints from ASE as a typical example
    ase_constraints = pyimport("ase.constraints")
    ase_io = pyimport("ase.io")
    structure = ase_io.read("artifacts/desorption_test.xyz", index=1)
    structure.set_constraint(ase_constraints.FixAtoms(indices=collect(0:17)))
    nqcd_atoms, nqcd_positions, nqcd_cell = NQCDynamics.convert_from_ase_atoms(structure)
    # Test the frozen model
    frozen_model = AdiabaticASEModel(structure)
    frozen_sim=Simulation(nqcd_atoms, frozen_model; cell=cell, temperature=10.0)

    # Pass: Atoms 1-18 are frozen, 19-56 are mobile
    @test NQCModels.mobileatoms(frozen_sim) == collect(19:56)
    test_dist = NQCDistributions.VelocityBoltzmann(10, frozen_sim)
    # Pass: First DOF of first atom should be frozen, so Dirac distribution
    findfirst(x->isa(x, Distributions.Dirac), test_dist.sampleable) == CartesianIndex(1,1)
    # Pass: First DOF of Atom 19 should be mobile, so Normal distribution
    findfirst(x->isa(x, Distributions.Normal), test_dist.sampleable) == CartesianIndex(1,19)
end
