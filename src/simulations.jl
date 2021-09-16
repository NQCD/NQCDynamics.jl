
using Unitful: @u_str
using UnitfulAtomic: austrip, auconvert

using .Calculators: AbstractCalculator, Calculator
using NonadiabaticModels: Model

abstract type AbstractSimulation{M,Calc<:AbstractCalculator,A<:Atoms,T,C<:AbstractCell} end

Base.broadcastable(sim::AbstractSimulation) = Ref(sim)

struct Simulation{M,Calc,A,T,C} <: AbstractSimulation{M,Calc,A,T,C}
    DoFs::Int
    temperature::T
    cell::C
    atoms::A
    calculator::Calc
    method::M
end

function Simulation(atoms::Atoms{S,T}, model::Model, method::M;
        DoFs::Integer=3, temperature=0u"K", cell::AbstractCell=InfiniteCell()) where {M,S,T}
    calc = Calculator(model, DoFs, length(atoms), T)
    Simulation(DoFs, temperature, cell, atoms, calc, method)
end

struct RingPolymerSimulation{M,Calc,A,T,C,B} <: AbstractSimulation{M,Calc,A,T,C}
    DoFs::Int
    temperature::T
    cell::C
    atoms::A
    calculator::Calc
    method::M
    beads::B
end

function RingPolymerSimulation(DoFs::Integer, temperature, cell::AbstractCell,
        atoms::Atoms{S,T}, model::Model, method::M,
        n_beads::Integer, quantum_nuclei::Vector{Symbol}) where {M,S,T}
        
    if isempty(quantum_nuclei)
        beads = RingPolymerParameters{T}(n_beads, get_temperature(temperature), length(atoms))
    else
        beads = RingPolymerParameters{T}(n_beads, get_temperature(temperature), atoms.types, quantum_nuclei)
    end
    
    calc = Calculator(model, DoFs, length(atoms), n_beads, T)
    RingPolymerSimulation(DoFs, temperature, cell, atoms, calc, method, beads)
end

function RingPolymerSimulation(atoms::Atoms, model::Model, method::M, n_beads::Integer;
        DoFs::Integer=3, temperature=0,
        cell::AbstractCell=InfiniteCell(), quantum_nuclei::Vector{Symbol}=Symbol[]) where {M}
    RingPolymerSimulation(DoFs, temperature, cell, atoms, model, method, n_beads, quantum_nuclei)
end

nbeads(sim::RingPolymerSimulation) = nbeads(sim.beads)

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
