module Dynamics

using ..Electronics
using ..Systems

using DiffEqBase

include("DynamicsVariables/DynamicsVariables.jl")
using .DynamicsVariables

export Phasespace
export MappingPhasespace
export SurfaceHoppingPhasespace

export get_momenta
export get_positions
export get_mapping_momenta
export get_mapping_positions
export get_density_matrix
export get_adiabatic_state

function differential!(du::DynamicalVariables, u::DynamicalVariables, p::AbstractSystem, t)
    get_positions(du) .= get_velocity(p, u)
    get_momenta(du) .= get_force(p, u)
end

function get_velocity(p::System, u::DynamicalVariables)
    get_momenta(u) ./ masses(p)
end

function get_velocity(p::RingPolymerSystem, u::DynamicalVariables)
    get_momenta(u) ./ get_bead_masses(n_beads(p), masses(p))
end

function get_force(p::System, u::DynamicalVariables)
    Electronics.calculate_derivative!(p.model, p.electronics, get_positions(u))
    -p.electronics.D0
end

function get_force(p::RingPolymerSystem, u::DynamicalVariables)
    Electronics.calculate_derivative!.(
        Ref(p.model),
        p.electronics,
        [get_positions(u)[i:n_beads(p):end] for i=1:n_beads(p)]
        )
    ground_state_force = -hcat([e.D0 for e in p.electronics]...)'[:]
    for (F, R, M) in zip(
        bead_iterator(n_beads(p), ground_state_force),
        get_bead_positions(u, n_beads(p)),
        get_bead_masses(n_beads(p), masses(p))
        )
        F .-= 2p.ring_polymer.springs * R * M
    end
    ground_state_force
end

end # module