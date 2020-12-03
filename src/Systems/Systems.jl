module Systems

using PeriodicTable
using Unitful
using UnitfulAtomic

using ..Atoms
using ..Models
using ..Electronics

include("ring_polymer.jl")

export AbstractSystem
export System
export RingPolymerSystem

export n_beads
export n_atoms
export masses

abstract type AbstractSystem end

"""
    System

Top level container for all the parametric quantities needed.
"""
struct System <: AbstractSystem
    n_DoF::UInt8
    atomic_parameters::AtomicParameters
    model::Models.Model
    electronics::Electronics.ElectronicContainer
end

function System(atomic_parameters::AtomicParameters, model::Models.Model, n_DoF::Integer=3)
    System(n_DoF, atomic_parameters, model, Electronics.ElectronicContainer(model.n_states, n_DoF*atomic_parameters.n_atoms))
end

struct RingPolymerSystem <: AbstractSystem
    atomic_parameters::AtomicParameters
    model::Models.Model
    electronics::Vector{Electronics.ElectronicContainer}
    ring_polymer::RingPolymerParameters
end

function RingPolymerSystem(atomic_parameters::AtomicParameters{T}, model::Models.Model, n_beads::Integer, temperature::Real, n_DoF::Integer=3) where {T<:AbstractFloat}
    electronics = [Electronics.ElectronicContainer{T}(model.n_states, n_DoF*atomic_parameters.n_atoms) for _=1:n_beads]
    ring_polymer = Systems.RingPolymerParameters{T}(n_beads, temperature)
    RingPolymerSystem(atomic_parameters, model, electronics, ring_polymer)
end

n_atoms(system::AbstractSystem) = system.atomic_parameters.n_atoms
masses(system::AbstractSystem) = repeat(system.atomic_parameters.masses, inner=system.n_DoF)
n_beads(system::RingPolymerSystem) = system.ring_polymer.n_beads

end # module


