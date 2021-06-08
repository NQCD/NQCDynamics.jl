using .Calculators: AbstractCalculator, Calculator
using Unitful
using UnitfulAtomic

export AbstractSimulation
export Simulation
export RingPolymerSimulation
export Method
export get_temperature

abstract type AbstractSimulation{M,Calc<:AbstractCalculator,A<:Atoms,T,C<:AbstractCell} end

Base.broadcastable(sim::AbstractSimulation) = Ref(sim)

struct Simulation{M,Calc,A,T,C} <: AbstractSimulation{M,Calc,A,T,C}
    DoFs::UInt8
    temperature::T
    cell::C
    atoms::A
    calculator::Calc
    method::M
end

function Simulation(DoFs::Integer, temperature, cell::AbstractCell,
        atoms::Atoms{S,T}, model::Model, method::M) where {M,S,T}
    calc = Calculator(model, DoFs, length(atoms), T)
    Simulation(UInt8(DoFs), temperature, cell, atoms, calc, method)
end

function Simulation(atoms::Atoms, model::Model, method::M;
        DoFs::Integer=3, temperature=0u"K",
        cell::AbstractCell=InfiniteCell()) where {M}
    Simulation(DoFs, temperature, cell, atoms, model, method)
end

function get_temperature(sim::AbstractSimulation, t::Real=0)
    t = auconvert(u"fs", t)
    get_temperature(sim.temperature, t)
end
get_temperature(temperature::Number, t=0) = austrip(temperature)
get_temperature(temperature::Function, t=0) = austrip(temperature(t))

struct RingPolymerSimulation{M,Calc,A,T,C,B} <: AbstractSimulation{M,Calc,A,T,C}
    DoFs::UInt8
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
    RingPolymerSimulation(UInt8(DoFs), temperature, cell, atoms, calc, method, beads)
end

function RingPolymerSimulation(atoms::Atoms, model::Model, method::M, n_beads::Integer;
        DoFs::Integer=3, temperature=0,
        cell::AbstractCell=InfiniteCell(), quantum_nuclei::Vector{Symbol}=Symbol[]) where {M}
    RingPolymerSimulation(DoFs, temperature, cell, atoms, model, method, n_beads, quantum_nuclei)
end
