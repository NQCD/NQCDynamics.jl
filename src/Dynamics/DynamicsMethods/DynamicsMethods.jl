module DynamicsMethods

using ..DynamicsVariables
using ....Atoms
using ....Systems
using ....Electronics
using ....Models

using Unitful
using UnitfulAtomic

export differential!

function differential!(du::DynamicalVariables, u::DynamicalVariables, p::AbstractSystem, t)
    set_velocity!(du, u, p)
    set_force!(du, u, p)
end

function set_velocity!(du::DynamicalVariables, u::DynamicalVariables, p::System{T}) where {T}
    get_positions(du) .= get_momenta(u) ./ masses(p)'
end

function set_velocity!(du::DynamicalVariables, u::DynamicalVariables, p::RingPolymerSystem{T}) where{T}
    for i=1:n_beads(p)
        get_positions(du, i) .= get_momenta(u, i) ./ masses(p)'
    end
end

function set_force!(du::DynamicalVariables, u::DynamicalVariables, p::System{T}) where {T}
    Electronics.calculate_derivative!(p.model, p.electronics, get_positions(u))
    get_momenta(du) .= -p.electronics.D0
end

function set_force!(du::DynamicalVariables, u::DynamicalVariables, p::RingPolymerSystem{T}) where {T}
    for i=1:n_beads(p)
        Electronics.calculate_derivative!(p.model, p.electronics[i], get_positions(u, i))
        get_momenta(du, i) .= -p.electronics[i].D0
    end
    apply_interbead_coupling!(du, u, p)
end

function apply_interbead_coupling!(du::DynamicalVariables, u::DynamicalVariables, p::RingPolymerSystem{T}) where {T}
    for i in p.ring_polymer.quantum_atoms
        for j=1:n_DoF(p)
            get_momenta(du)[j,i,:] .-= 2p.ring_polymer.springs*get_positions(u)[j,i,:] .* masses(p)[i]
        end
    end
end

include("langevin.jl")
export Langevin
include("mdef.jl")
export MDEF

include("IESH_dynamics.jl")
export IESH_dynamics

end # module