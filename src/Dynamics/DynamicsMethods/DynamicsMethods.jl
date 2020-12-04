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
    get_positions(du) .= get_momenta(u) ./ masses(p)
end

function set_velocity!(du::DynamicalVariables, u::DynamicalVariables, p::RingPolymerSystem{T}) where{T}
    get_positions(du) .= hcat([P ./ masses(p) for P in eachcol(get_momenta(u))]...)
end

function set_force!(du::DynamicalVariables, u::DynamicalVariables, p::System{T}) where {T}
    Electronics.calculate_derivative!(p.model, p.electronics, get_positions(u))
    get_momenta(du) .= -p.electronics.D0
end

function set_force!(du::DynamicalVariables, u::DynamicalVariables, p::RingPolymerSystem{T}) where {T}
    Electronics.calculate_derivative!.(Ref(p.model), p.electronics, eachcol(get_positions(u)))

    get_momenta(du) .= -hcat([e.D0 for e in p.electronics]...)
    apply_interbead_coupling!(du, u, p)
end

function apply_interbead_coupling!(du::DynamicalVariables, u::DynamicalVariables, p::RingPolymerSystem{T}) where {T}
    for (i, bead) in enumerate(eachrow(get_positions(u)))
        get_momenta(du)[i,:] .-= 2p.ring_polymer.springs * bead * masses(p)[i]
    end
end

include("langevin.jl")
export Langevin

end # module