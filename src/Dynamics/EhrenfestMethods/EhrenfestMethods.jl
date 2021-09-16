
module EhrenfestMethods

"""
Abstract type for Ehrenfest method.
"""
abstract type AbstractEhrenfest <: Method end

function motion!(du, u, sim::AbstractSimulation{<:AbstractEhrenfest}, t)
    dr = get_positions(du)
    dv = get_velocities(du)
    dσ = get_quantum_subsystem(du)

    r = get_positions(u)
    v = get_velocities(u)
    σ = get_quantum_subsystem(u)

    velocity!(dr, v, r, sim, t)
    Calculators.update_electronics!(sim.calculator, r)
    acceleration!(dv, v, r, sim, t, σ)
    set_quantum_derivative!(dσ, v, σ, sim)
end

include("ehrenfest.jl")
include("ehrenfest_rpmd.jl")

end # module
