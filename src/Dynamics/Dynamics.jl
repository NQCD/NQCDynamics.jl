module Dynamics

using ..Electronics
using ..Systems

using DiffEqBase

include("DynamicsVariables/DynamicsVariables.jl")
using .DynamicsVariables

export Phasespace
export MappingPhasespace
export SurfaceHoppingPhasespace
export RingPolymerPhasespace

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

get_velocity(p::System, u::DynamicalVariables) = get_momenta(u) ./ masses(p)

function get_velocity(p::RingPolymerSystem, u::RingPolymerDynamicalVariables)
    hcat([P ./ masses(p) for P in eachcol(get_momenta(u))]...)
end

function get_force(p::System, u::DynamicalVariables)
    Electronics.calculate_derivative!(p.model, p.electronics, get_positions(u))
    -p.electronics.D0
end

function get_force(p::RingPolymerSystem, u::RingPolymerDynamicalVariables)
    Electronics.calculate_derivative!.(Ref(p.model), p.electronics, eachcol(get_positions(u)))

    ground_state_force = -hcat([e.D0 for e in p.electronics]...)
    for (i, bead) in enumerate(eachrow(get_positions(u)))
        ground_state_force[i,:] .-= 2p.ring_polymer.springs * bead * masses(p)[i]
    end

    ground_state_force
end

end # module