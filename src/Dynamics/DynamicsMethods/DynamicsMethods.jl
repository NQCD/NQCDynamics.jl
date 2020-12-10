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
    dofs = n_DoF(p)*n_atoms(p)
    atom_index = 1
    for i=1:dofs
        # Multiplication with the view of positions doesn't work due to a bug with the ArrayPartition
        get_bead_momenta(du, i, dofs) .-= 2p.ring_polymer.springs * get_positions(u)[i:dofs:end] .* masses(p)[atom_index]
        if i % n_DoF(p) == 0; atom_index += 1 end
    end
end

include("langevin.jl")
export Langevin

end # module