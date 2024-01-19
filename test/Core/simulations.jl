using Test
using NQCDynamics
using Unitful, UnitfulAtomic

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

@testset "Thermostats" begin
    # Init thermostat with unitful quality
    thermostat1=Thermostat(10u"K", [1,2])
    @test NQCDynamics.get_temperature(thermostat1, 0) == NQCDynamics.get_temperature(thermostat1, 100)
    # Init thermostat with function
    thermostat2=Thermostat(x->(10+x*u"fs^-1")*u"K", [2,3])
    @test NQCDynamics.get_temperature(thermostat2) == austrip(10u"K")
    @test NQCDynamics.get_temperature(thermostat2, 10u"fs") == austrip(20u"K")
    # Test that thermostats with overlapping indices throw an error
    @test_throws DomainError Simulation(atoms, model; temperature=[thermostat1, thermostat2])
    # Test system size / total thermostat size mismatch
    @test_throws DomainError Simulation(Atoms([:N, :H, :H, :H,]), model; temperature=thermostat1)
    thermostat3=thermostat2=Thermostat(x->(10+x*u"fs^-1")*u"K", [3,4])
    # Combined temperature evaluation
    sim = Simulation(Atoms([:N, :H, :H, :H,]), model; temperature=[thermostat1, thermostat3])
    @test NQCDynamics.get_temperature(sim, austrip(1u"fs")) == austrip.([10u"K", 10u"K", 11u"K", 11u"K"])
end

@testset "get_temperature(Simulation)" begin
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
