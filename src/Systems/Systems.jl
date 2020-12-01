module Systems

using PeriodicTable
using Unitful
using UnitfulAtomic

using ..Models
using ..Electronics

include("cell.jl")
include("atomic_parameters.jl")

export System

abstract type AbstractSystem end

"""
    System

Top level container for all the parametric quantities needed.
"""
struct System <: AbstractSystem
    atomic_parameters::AtomicParameters
    model::Models.Model
    electronics::Electronics.ElectronicContainer
end

function System(atomic_parameters::AtomicParameters, model::Models.Model)
    System(atomic_parameters, model, Electronics.ElectronicContainer(model.n_states, atomic_parameters.n_atoms))
end

struct RingPolymerSystem <: AbstractSystem
    atomic_parameters::AtomicParameters
    model::Models.Model
    electronics::Electronics.ElectronicContainer
end

end # module


