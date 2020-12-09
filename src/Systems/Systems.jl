module Systems

using PeriodicTable
using Unitful
using UnitfulAtomic

using ..Atoms
using ..Models
using ..Electronics

export AbstractSystem
export System
export RingPolymerSystem
export Classical

export n_beads
export n_atoms
export masses
export n_DoF

abstract type DynamicsParameters end
abstract type AbstractSystem{D<:DynamicsParameters} end

struct Classical <: DynamicsParameters end

include("ring_polymer.jl")

"""
    System

Top level container for all the parametric quantities needed.
"""
struct System{D} <: AbstractSystem{D}
    n_DoF::UInt8
    atomic_parameters::AtomicParameters
    model::Models.Model
    electronics::Electronics.ElectronicContainer
    dynamics::D
    function System{D}(n_DoF::Integer, atomic_parameters::AtomicParameters, model::Models.Model,
        dynamics::D) where {D<:DynamicsParameters}
        electronics = Electronics.ElectronicContainer(model.n_states, n_DoF, atomic_parameters.n_atoms)
        new{D}(n_DoF, atomic_parameters, model, electronics, dynamics)
    end
end

function System(atomic_parameters::AtomicParameters, model::Models.Model, n_DoF::Integer=3)
    System{Classical}(n_DoF, atomic_parameters, model, Classical())
end

struct RingPolymerSystem{D} <: AbstractSystem{D}
    n_DoF::UInt8
    atomic_parameters::AtomicParameters
    model::Models.Model
    electronics::Vector{Electronics.ElectronicContainer}
    ring_polymer::RingPolymerParameters
    dynamics::D
end

function RingPolymerSystem(atomic_parameters::AtomicParameters{T}, model::Models.Model, n_beads::Integer, temperature::Real, n_DoF::Integer=3) where {T<:AbstractFloat}
    electronics = [Electronics.ElectronicContainer{T}(model.n_states, n_DoF, atomic_parameters.n_atoms) for _=1:n_beads]
    ring_polymer = Systems.RingPolymerParameters{T}(n_beads, temperature)
    RingPolymerSystem{Classical}(n_DoF, atomic_parameters, model, electronics, ring_polymer, Classical())
end

n_atoms(system::AbstractSystem) = system.atomic_parameters.n_atoms
masses(system::AbstractSystem) = system.atomic_parameters.masses
n_beads(system::RingPolymerSystem) = system.ring_polymer.n_beads
n_DoF(system::RingPolymerSystem) = system.n_DoF

end # module


