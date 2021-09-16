using Test
using NonadiabaticMolecularDynamics
using Unitful

atoms = Atoms([:C, :H])
cell = InfiniteCell()

model = NonadiabaticModels.Free()

@test RingPolymerSimulation(3, 100, cell, atoms, model, Classical(), 10, [:H]) isa RingPolymerSimulation

@test Simulation(atoms, model, Classical()) isa Simulation
@test RingPolymerSimulation(atoms, model, Classical(), 10) isa RingPolymerSimulation

@test Simulation(atoms, model) isa Simulation{Classical}
@test Simulation{Classical}(atoms, model) isa Simulation{Classical}
@test RingPolymerSimulation{Classical}(atoms, model, 10) isa RingPolymerSimulation{Classical}

# @test Simulation{MDEF}(atoms, model) isa Simulation{<:MDEF}

@test RingPolymerSimulation{NRPMD}(atoms, NonadiabaticModels.DoubleWell(), 10) isa RingPolymerSimulation{<:NRPMD}
@test Simulation{Langevin}(atoms, model) isa Simulation{<:Langevin}
@test Simulation{FSSH}(atoms, NonadiabaticModels.DoubleWell()) isa Simulation{<:FSSH}

@testset "get_temperature" begin
    sim = Simulation(atoms, model; temperature=10u"K")
    @test NonadiabaticMolecularDynamics.get_temperature(sim) isa Real
    sim = Simulation(atoms, model; temperature=x->10u"K")
    @test NonadiabaticMolecularDynamics.get_temperature(sim) isa Real
    sim = Simulation(atoms, model; temperature=10)
    @test NonadiabaticMolecularDynamics.get_temperature(sim, 1) isa Real
end
