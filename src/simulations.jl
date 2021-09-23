
using Unitful: @u_str
using UnitfulAtomic: austrip, auconvert

using .Calculators: AbstractCalculator, Calculator
using NonadiabaticModels: Model

abstract type AbstractSimulation{M,Calc<:AbstractCalculator,A<:Atoms,T,C<:AbstractCell} end

Base.broadcastable(sim::AbstractSimulation) = Ref(sim)

struct Simulation{M,Calc,A,T,C} <: AbstractSimulation{M,Calc,A,T,C}
    temperature::T
    cell::C
    atoms::A
    calculator::Calc
    method::M
end

function Simulation(atoms::Atoms{S,T}, model::Model, method::M;
        temperature=0u"K", cell::AbstractCell=InfiniteCell()) where {M,S,T}
    calc = Calculator(model, length(atoms), T)
    Simulation(temperature, cell, atoms, calc, method)
end

struct RingPolymerSimulation{M,Calc,A,T,C,B} <: AbstractSimulation{M,Calc,A,T,C}
    temperature::T
    cell::C
    atoms::A
    calculator::Calc
    method::M
    beads::B
end

function RingPolymerSimulation(temperature, cell::AbstractCell,
        atoms::Atoms{S,T}, model::Model, method::M,
        n_beads::Integer, quantum_nuclei::Vector{Symbol}) where {M,S,T}
        
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

NonadiabaticModels.nstates(sim::AbstractSimulation) = NonadiabaticModels.nstates(sim.calculator.model)

NonadiabaticModels.ndofs(sim::AbstractSimulation) = NonadiabaticModels.ndofs(sim.calculator.model)
natoms(sim::AbstractSimulation) = length(sim.atoms)
RingPolymers.nbeads(sim::RingPolymerSimulation) = RingPolymers.nbeads(sim.beads)
masses(sim::AbstractSimulation) = sim.atoms.masses
masses(sim::AbstractSimulation, i) = sim.atoms.masses[i]

Base.size(sim::Simulation) = (ndofs(sim), natoms(sim))
Base.size(sim::RingPolymerSimulation) = (ndofs(sim), natoms(sim), RingPolymers.nbeads(sim))

function get_temperature(sim::Simulation, t::Real=0)
    t = auconvert(u"fs", t)
    get_temperature(sim.temperature, t)
end
function get_temperature(sim::RingPolymerSimulation, t::Real=0)
    t = auconvert(u"fs", t)
    get_temperature(sim.temperature, t) * length(sim.beads)
end
get_temperature(temperature::Number, t=0) = austrip(temperature)
get_temperature(temperature::Function, t=0) = austrip(temperature(t))

function Base.show(io::IO, sim::Simulation{M}) where {M}
    print(io, "Simulation{$M}:\n  ", sim.atoms, "\n  ", sim.calculator.model)
end

function Base.show(io::IO, sim::RingPolymerSimulation{M}) where {M}
    print(io, "RingPolymerSimulation{$M}:\n\n  ", sim.atoms, "\n\n  ", sim.calculator.model,
          "\n  with ", length(sim.beads), " beads.")
end
