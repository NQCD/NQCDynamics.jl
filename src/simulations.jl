using .Calculators: AbstractCalculator, Calculator
using .Models: Model
using Unitful
using UnitfulAtomic

export AbstractSimulation
export Simulation
export RingPolymerSimulation
export Method

abstract type AbstractSimulation{M,S,T<:AbstractFloat,C<:AbstractCell,A<:AbstractCalculator} end

struct Simulation{M,S,T,C,A} <: AbstractSimulation{M,S,T,C,A}
    DoFs::UInt8
    temperature::T
    cell::C
    atoms::Atoms{S,T}
    calculator::A
    method::M
    function Simulation(DoFs::Integer, temperature::Real, cell::AbstractCell{T},
            atoms::Atoms{S,T}, model::Model, method::M) where {M,S,T}
        calc = Calculator(model, DoFs, length(atoms), T)
        new{M,S,T,typeof(cell),typeof(calc)}(DoFs, temperature, cell, atoms, calc, method)
    end
end

function Simulation(DoFs::Integer, temperature::Unitful.Temperature, cell::AbstractCell{T},
        atoms::Atoms{S,T}, model::Model, method::M) where {M,S,T}
    Simulation(DoFs, austrip(temperature), cell, atoms, model, method)
end

function Simulation(atoms::Atoms{S,T}, model::Model, method::M;
        DoFs::Integer=3, temperature::Unitful.Temperature=0u"K",
        cell::AbstractCell{T}=InfiniteCell{T}()) where {M,S,T}
    Simulation(DoFs, temperature, cell, atoms, model, method)
end

struct RingPolymerSimulation{M,S,T,C,A} <: AbstractSimulation{M,S,T,C,A}
    DoFs::UInt8
    temperature::T
    cell::C
    atoms::Atoms{S,T}
    calculator::A
    method::M
    beads::RingPolymerParameters{T}
    function RingPolymerSimulation(DoFs::Integer, temperature::Real, cell::AbstractCell{T},
            atoms::Atoms{S,T}, model::Model, method::M,
            n_beads::Integer, quantum_nuclei::Vector{Symbol}=Symbol[]) where {M,S,T}
            
        if isempty(quantum_nuclei)
            beads = RingPolymerParameters{T}(n_beads, temperature, length(atoms))
        else
            beads = RingPolymerParameters{T}(n_beads, temperature, atoms.types, quantum_nuclei)
        end
        
        calc = Calculator(model, DoFs, length(atoms), n_beads, T)
        new{M,S,T,typeof(cell),typeof(calc)}(DoFs, temperature, cell, atoms, calc, method, beads)
    end
end

function RingPolymerSimulation(DoFs::Integer, temperature::Unitful.Temperature, cell::AbstractCell{T},
        atoms::Atoms{S,T}, model::Model, method::M,
        n_beads::Integer, quantum_nuclei::Vector{Symbol}=Symbol[]) where {M,S,T}
    RingPolymerSimulation(DoFs, austrip(temperature), cell, atoms, model, method, n_beads, quantum_nuclei)
end

function RingPolymerSimulation(atoms::Atoms{S,T}, model::Model, method::M, n_beads::Integer;
        DoFs::Integer=3, temperature::Unitful.Temperature=0u"K",
        cell::AbstractCell{T}=InfiniteCell{T}(), quantum_nuclei::Vector{Symbol}=Symbol[]) where {M,S,T}
    RingPolymerSimulation(DoFs, austrip(temperature), cell, atoms, model, method, n_beads, quantum_nuclei)
end