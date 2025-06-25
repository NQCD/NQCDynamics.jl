
using Unitful: @u_str
using UnitfulAtomic: austrip, auconvert

using NQCCalculators
using NQCModels: Model, Subsystem, CompositeModel

abstract type AbstractSimulation{M} end

Base.broadcastable(sim::AbstractSimulation) = Ref(sim)

struct Simulation{M,Cache,T,Ttype,C,S} <: AbstractSimulation{M}
    temperature::Ttype
    cell::C
    atoms::Atoms{T}
    cache::Cache
    method::M
    solver::S
end

"""
    Simulation(atoms::Atoms{T}, model::Model, method::M;
        temperature=0u"K", cell::AbstractCell=InfiniteCell()) where {M,S,T}

Simulation parameters that controls the types of atoms, interactions,
dynamics method, temperature and simulation cell.
"""
function Simulation(atoms::Atoms{T}, model::Model, method::M;
        temperature=0u"K", cell::AbstractCell=InfiniteCell(), solver::Symbol=:exact) where {M,T}
    cache = Create_Cache(model, length(atoms), T)
    # If a thermostat is provided, check it covers the whole system. 
    isa(temperature, TemperatureSetting{Vector{Int}}) ? throw(DomainError(temperature, "TemperatureSetting must apply to all atoms.")) : nothing
    # If multiple TemperatureSettings are provided, check that each atom only has one thermostat applied to it.
    if isa(temperature, Vector{<:TemperatureSetting})
        indices = vcat([thermostat.indices for thermostat in temperature]...)
        if length(unique(indices)) != length(atoms.masses)
            throw(DomainError(temperature, "Every atom must have a TemperatureSetting applied to it."))
        end
        if length(indices) != length(unique(indices))
            throw(DomainError(temperature, "Atoms can only have one thermostat applied to them."))
        end
    end
    Simulation(temperature, cell, atoms, cache, method, solver)
end

struct RingPolymerSimulation{M,Cache,T,Ttype,C,B,S} <: AbstractSimulation{M}
    temperature::Ttype
    cell::C
    atoms::Atoms{T}
    cache::Cache
    method::M
    beads::B
    solver::S
end

function RingPolymerSimulation(temperature, cell::AbstractCell,
        atoms::Atoms{T}, model::Model, method::M,
        n_beads::Integer, quantum_nuclei::Vector{Symbol}; solver::Symbol=:exact) where {M,T}
        
    if isempty(quantum_nuclei)
        beads = RingPolymers.RingPolymerParameters{T}(n_beads, get_temperature(temperature), length(atoms))
    else
        beads = RingPolymers.RingPolymerParameters{T}(n_beads, get_temperature(temperature), atoms.types, quantum_nuclei)
    end
    
    cache = Create_Cache(model, length(atoms), n_beads, T)
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
get_temperature(temperature::Function, t=0) = austrip(temperature(t))

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

"""
A TemperatureSetting contains both temperature information and the atom indices within a Simulation that it is applied to. 

**If you don't need to apply different temperatures to different parts of your Simulation, assign the `temperature` keyword argument of your simulation to a number or function.**

## Parameters

`value`: A temperature function. This can be a Number type for constant temperatures, or a function taking the time in Unitful `u"ps"` as input and giving a temperature in Unitful `u"K"` as output. 

`indices`: Indices of the atoms to which this thermostat is applied. Can be a range of indices, a single `Int`, or a `Vector{Int}`.

"""
struct TemperatureSetting{I}
    value::Function
    indices::I
end

function Base.show(io::IO, temperature::TemperatureSetting)
    print(io, "TemperatureSetting:\n\tT(t=0) = $(get_temperature(temperature, 0u"fs"))\n\tApplies to atoms: $(temperature.indices)\n")
end

function TemperatureSetting(value, indices=:)
    if isa(indices, UnitRange)
        indices=collect(indices)
    elseif isa(indices, Int)
        indices=[indices]
    end
    if !isa(value, Function)
        temperature_function(t)=value
        return TemperatureSetting(temperature_function, indices)
    else
        return TemperatureSetting(value, indices)
    end
end

"""
    TemperatureSetting(temperature, subsystem::NQCModels.Subsystem)

Apply a `TemperatureSetting` to all atoms in a `Subsystem`. 
"""
function TemperatureSetting(temperature, subsystem::NQCModels.Subsystem)
    TemperatureSetting(temperature, subsystem.indices)
end

get_temperature(thermostat::TemperatureSetting{<:Any}, t=0u"fs") = austrip(thermostat.value(t))

"""
    get_temperature(thermostats::Vector{<:TemperatureSetting}, t=0u"fs")

Gets the temperature from multiple `TemperatureSetting`s and returns a vector of the temperature applied to each atom. 
"""
function get_temperature(thermostats::Vector{<:TemperatureSetting}, t=0u"fs")
    indices=vcat([thermostat.indices for thermostat in thermostats]...)
    temperature_vector=zeros(Number, length(indices))
    for thermostat in thermostats
        temperature_vector[thermostat.indices] .= austrip(thermostat.value(t))
    end
    return temperature_vector
end

