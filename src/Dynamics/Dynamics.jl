module Dynamics

using ..Electronics
using ..Systems

using RecursiveArrayTools
using DiffEqBase
using Unitful
using UnitfulAtomic
include("DynamicalVariables/phasespace.jl")

function differential!(du::DynamicalVariables, u::DynamicalVariables, p::System, t)
    get_positions(du) .= get_momenta(u) ./ p.atomic_parameters.masses
    get_momenta(du) .= get_force(p, u)
end

function get_force(p::System, u::DynamicalVariables)
    Electronics.calculate_derivative!(p.model, p.electronics, get_positions(u))
    -p.electronics.D0
end

end # module