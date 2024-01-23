
using Unitful: @u_str
using UnitfulAtomic: austrip, auconvert

using .Calculators: AbstractCalculator, Calculator
using NQCModels: Model, Subsystem, CompositeModel

abstract type AbstractSimulation{M} end

Base.broadcastable(sim::AbstractSimulation) = Ref(sim)

struct Simulation{M,Calc,T,Ttype,C} <: AbstractSimulation{M}
    temperature::Ttype
    cell::C
    atoms::Atoms{T}
    calculator::Calc
    method::M
end

"""
    Simulation(atoms::Atoms{T}, model::Model, method::M;
        temperature=0u"K", cell::AbstractCell=InfiniteCell()) where {M,S,T}

Simulation parameters that controls the types of atoms, interactions,
dynamics method, temperature and simulation cell.
"""
function Simulation(atoms::Atoms{T}, model::Model, method::M;
        temperature=0u"K", cell::AbstractCell=InfiniteCell()) where {M,T}
    calc = Calculator(model, length(atoms), T)
    # If a thermostat is provided, check it covers the whole system. 
    isa(temperature, Thermostat{Vector{Int}}) ? throw(DomainError(temperature, "Thermostat must apply to all atoms.")) : nothing
    # If multiple Thermostats are provided, check that each atom only has one thermostat applied to it.
    if isa(temperature, Vector{<:Thermostat})
        indices = vcat([thermostat.indices for thermostat in temperature]...)
        if length(unique(indices)) != length(atoms.masses)
            throw(DomainError(temperature, "Every atom must have a Thermostat applied to it."))
        end
        if length(indices) != length(unique(indices))
            throw(DomainError(temperature, "Atoms can only have one thermostat applied to them."))
        end
    end
    Simulation(temperature, cell, atoms, calc, method)
end

struct RingPolymerSimulation{M,Calc,T,Ttype,C,B} <: AbstractSimulation{M}
    temperature::Ttype
    cell::C
    atoms::Atoms{T}
    calculator::Calc
    method::M
    beads::B
end

function RingPolymerSimulation(temperature, cell::AbstractCell,
        atoms::Atoms{T}, model::Model, method::M,
        n_beads::Integer, quantum_nuclei::Vector{Symbol}) where {M,T}
        
    if isempty(quantum_nuclei)
        beads = RingPolymers.RingPolymerParameters{T}(n_beads, get_temperature(temperature), length(atoms))
    else
        beads = RingPolymers.RingPolymerParameters{T}(n_beads, get_temperature(temperature), atoms.types, quantum_nuclei)
    end
    
    calc = Calculator(model, length(atoms), n_beads, T)
    RingPolymerSimulation(temperature, cell, atoms, calc, method, beads)
end

function RingPolymerSimulation(atoms::Atoms, model::Model, method::M, n_beads::Integer;
        temperature=0,
        cell::AbstractCell=InfiniteCell(), quantum_nuclei::Vector{Symbol}=Symbol[]) where {M}
    RingPolymerSimulation(temperature, cell, atoms, model, method, n_beads, quantum_nuclei)
end

NQCModels.nstates(sim::AbstractSimulation) = NQCModels.nstates(sim.calculator)
NQCModels.eachstate(sim::AbstractSimulation) = NQCModels.eachstate(sim.calculator)
NQCModels.nelectrons(sim::AbstractSimulation) = NQCModels.nelectrons(sim.calculator)
NQCModels.eachelectron(sim::AbstractSimulation) = NQCModels.eachelectron(sim.calculator)
NQCModels.mobileatoms(sim::AbstractSimulation) = NQCModels.mobileatoms(sim.calculator)
NQCModels.dofs(sim::AbstractSimulation) = NQCModels.dofs(sim.calculator)
NQCModels.fermilevel(sim::AbstractSimulation) = NQCModels.fermilevel(sim.calculator)

NQCModels.ndofs(sim::AbstractSimulation) = NQCModels.ndofs(sim.calculator.model)
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
    print(io, "Simulation{$M}:\n  ", sim.atoms, "\n  ", sim.calculator.model)
end

function Base.show(io::IO, sim::RingPolymerSimulation{M}) where {M}
    print(io, "RingPolymerSimulation{$M}:\n\n  ", sim.atoms, "\n\n  ", sim.calculator.model,
          "\n  with ", length(sim.beads), " beads.")
end

"""
A Thermostat is defined by a temperature (either constant or time-dependent) as the atom indices within a structure that it is applied to. 

"""
struct Thermostat{I}
    temperature::Function
    indices::I
end

function Base.show(io::IO, thermostat::Thermostat)
    print(io, "Thermostat:\n\tT(0) = $(get_temperature(thermostat, 0u"fs"))\n\tApplies to atoms: $(thermostat.indices)\n")
end

function Thermostat(temperature, indices=:)
    if isa(indices, UnitRange)
        indices=collect(indices)
    elseif isa(indices, Int)
        indices=[indices]
    end
    if !isa(temperature, Function)
        temperature_function(t)=temperature
        return Thermostat(temperature_function, indices)
    else
        return Thermostat(temperature, indices)
    end
end

function Thermostat(temperature, subsystem::NQCModels.Subsystem)
    Thermostat(temperature, subsystem.indices)
end

get_temperature(thermostat::Thermostat{<:Any}, t=0u"fs") = austrip(thermostat.temperature(t))

function get_temperature(thermostats::Vector{<:Thermostat}, t=0u"fs")
    indices=vcat([thermostat.indices for thermostat in thermostats]...)
    temperature_vector=zeros(Number, length(indices))
    for thermostat in thermostats
        temperature_vector[thermostat.indices] .= austrip(thermostat.temperature(t))
    end
    return temperature_vector
end

