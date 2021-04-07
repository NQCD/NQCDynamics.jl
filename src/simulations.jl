using .Calculators: AbstractCalculator, Calculator
using .Models: Model
using Unitful
using UnitfulAtomic

export AbstractSimulation
export Simulation
export RingPolymerSimulation
export Method
export get_temperature

abstract type AbstractSimulation{M,Calc<:AbstractCalculator,A<:Atoms,T,C<:AbstractCell} end

struct Simulation{M,Calc,A,T,C} <: AbstractSimulation{M,Calc,A,T,C}
    DoFs::UInt8
    temperature::T
    cell::C
    atoms::A
    calculator::Calc
    method::M
end

Base.broadcastable(sim::AbstractSimulation) = Ref(sim)

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

function get_temperature(sim::AbstractSimulation{M,Calc,A,T,C}, t::Real=0) where {M,Calc,A,T,C}
    austrip(sim.temperature)
end
function get_temperature(sim::AbstractSimulation{M,Calc,A,T,C}, t::Real=0) where {M,Calc,A,T<:Function,C}
    t = auconvert(u"fs", t)
    austrip(sim.temperature(t))
end

struct RingPolymerSimulation{M,Calc,A,T,C,B} <: AbstractSimulation{M,Calc,A,T,C}
    DoFs::UInt8
    temperature::T
    cell::C
    atoms::A
    calculator::Calc
    method::M
    beads::B
    function RingPolymerSimulation(DoFs::Integer, temperature::Real, cell::AbstractCell,
            atoms::Atoms{S,T}, model::Model, method::M,
            n_beads::Integer, quantum_nuclei::Vector{Symbol}=Symbol[]) where {M,S,T}
            
        if isempty(quantum_nuclei)
            beads = RingPolymerParameters{T}(n_beads, temperature, length(atoms))
        else
            beads = RingPolymerParameters{T}(n_beads, temperature, atoms.types, quantum_nuclei)
        end
        
        calc = Calculator(model, DoFs, length(atoms), n_beads, T)
        new{M,typeof(calc),typeof(atoms),typeof(temperature),typeof(cell),typeof(beads)}(DoFs, temperature, cell, atoms, calc, method, beads)
    end
end

function RingPolymerSimulation(DoFs::Integer, temperature::Unitful.Temperature, cell::AbstractCell,
        atoms::Atoms, model::Model, method::M,
        n_beads::Integer, quantum_nuclei::Vector{Symbol}=Symbol[]) where {M}
    RingPolymerSimulation(DoFs, austrip(temperature), cell, atoms, model, method, n_beads, quantum_nuclei)
end

function RingPolymerSimulation(atoms::Atoms, model::Model, method::M, n_beads::Integer;
        DoFs::Integer=3, temperature::Unitful.Temperature=0u"K",
        cell::AbstractCell=InfiniteCell(), quantum_nuclei::Vector{Symbol}=Symbol[]) where {M}
    RingPolymerSimulation(DoFs, austrip(temperature), cell, atoms, model, method, n_beads, quantum_nuclei)
end
