
using Unitful: @u_str
using UnitfulAtomic: austrip, auconvert

using NQCCalculators
using NQCModels: Model, Subsystem, CompositeModel
import Adapt

abstract type AbstractSimulation{M} end

Base.broadcastable(sim::AbstractSimulation) = Ref(sim)

struct Simulation{M,Cache,A<:Atoms,Ttype,C,S} <: AbstractSimulation{M}
    temperature::Ttype
    cell::C
    atoms::A
    cache::Cache
    method::M
    solver::S
end

Adapt.@adapt_structure Simulation

"""
    Simulation(atoms::Atoms, model::Model, method::M;
        temperature=0u"K", cell::AbstractCell=InfiniteCell()) where {M,S}

Simulation parameters that controls the types of atoms, interactions,
dynamics method, temperature and simulation cell.
"""
function Simulation(atoms::Atoms, model::Model, method::M;
        temperature=0u"K", cell::AbstractCell=InfiniteCell(), solver::Symbol=:exact, device = NQCBase.CPUDevice) where {M}
    cache = Create_Cache(model, length(atoms), atoms.masses)
    Simulation(temperature, cell, atoms, cache, method, solver)
end

struct RingPolymerSimulation{M,Cache,A<:Atoms,Ttype,C,B,S} <: AbstractSimulation{M}
    temperature::Ttype
    cell::C
    atoms::A
    cache::Cache
    method::M
    beads::B
    solver::S
end

function RingPolymerSimulation(temperature, cell::AbstractCell,
        atoms::Atoms, model::Model, method::M,
        n_beads::Integer, quantum_nuclei::Vector{Symbol}; solver::Symbol=:exact) where {M}
    T = eltype(atoms.masses)

    if isempty(quantum_nuclei)
        beads = RingPolymers.RingPolymerParameters{T}(n_beads, get_temperature(temperature), length(atoms))
    else
        beads = RingPolymers.RingPolymerParameters{T}(n_beads, get_temperature(temperature), atoms.types, quantum_nuclei)
    end

    cache = Create_Cache(model, length(atoms), n_beads, atoms.masses)
    RingPolymerSimulation(temperature, cell, atoms, cache, method, beads, solver)
end

function RingPolymerSimulation(atoms::Atoms, model::Model, method::M, n_beads::Integer; temperature=0,
        cell::AbstractCell=InfiniteCell(), quantum_nuclei::Vector{Symbol}=Symbol[], solver::Symbol=:exact) where {M}
    RingPolymerSimulation(temperature, cell, atoms, model, method, n_beads, quantum_nuclei; solver=solver)
end

NQCModels.nstates(sim::AbstractSimulation) = NQCModels.nstates(sim.cache)
NQCModels.eachstate(sim::AbstractSimulation) = NQCModels.eachstate(sim.cache)
NQCModels.nelectrons(sim::AbstractSimulation) = NQCModels.nelectrons(sim.cache)
NQCModels.eachelectron(sim::AbstractSimulation) = NQCModels.eachelectron(sim.cache)
NQCModels.mobileatoms(sim::AbstractSimulation) = NQCModels.mobileatoms(sim.cache)
NQCModels.dofs(sim::AbstractSimulation) = NQCModels.dofs(sim.cache)
NQCModels.fermilevel(sim::AbstractSimulation) = NQCModels.fermilevel(sim.cache)

NQCModels.ndofs(sim::AbstractSimulation) = NQCModels.ndofs(sim.cache.model)
natoms(sim::AbstractSimulation) = length(sim.atoms)
RingPolymers.nbeads(sim::RingPolymerSimulation) = RingPolymers.nbeads(sim.beads)
masses(sim::AbstractSimulation) = sim.atoms.masses
masses(sim::AbstractSimulation, i::Int) = sim.atoms.masses[i]

function masses(sim::AbstractSimulation, I::CartesianIndex)
    masses(sim, I[2])
end

Base.size(sim::Simulation) = (ndofs(sim), natoms(sim))
Base.size(sim::RingPolymerSimulation) = (ndofs(sim), natoms(sim), RingPolymers.nbeads(sim))

function get_temperature(sim::AbstractSimulation, t::Real=0)
    t = auconvert(u"fs", t)
    get_temperature(sim.temperature, t)
end

get_temperature(temperature::Number, t=0) = austrip(temperature)
get_temperature(temperature::Function, t=0) = austrip.(temperature(t)) # return can be a <:Number or a length(Atoms) array.

function get_ring_polymer_temperature(sim::RingPolymerSimulation, t::Real=0)
    return @. get_temperature(sim, t) * nbeads(sim)
end

function Base.show(io::IO, sim::Simulation{M}) where {M}
    print(io, "Simulation{$M}:\n  ", sim.atoms, "\n  ", sim.cache.model)
end

function Base.show(io::IO, sim::RingPolymerSimulation{M}) where {M}
    print(io, "RingPolymerSimulation{$M}:\n\n  ", sim.atoms, "\n\n  ", sim.cache.model,
          "\n  with ", length(sim.beads), " beads.")
end
