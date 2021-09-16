
module EhrenfestMethods

using ....NonadiabaticMolecularDynamics:
    NonadiabaticMolecularDynamics,
    AbstractSimulation,
    Simulation,
    RingPolymerSimulation
using ....Calculators: Calculators
using ..Dynamics: Dynamics

"""
Abstract type for Ehrenfest method.
"""
abstract type AbstractEhrenfest <: Dynamics.Method end

function Dynamics.motion!(du, u, sim::AbstractSimulation{<:AbstractEhrenfest}, t)
    dr = get_positions(du)
    dv = get_velocities(du)
    dσ = get_quantum_subsystem(du)

    r = get_positions(u)
    v = get_velocities(u)
    σ = get_quantum_subsystem(u)

    DynamicsUtils.velocity!(dr, v)
    Calculators.update_electronics!(sim.calculator, r)
    acceleration!(dv, v, r, sim, t, σ)
    Dynamics.set_quantum_derivative!(dσ, v, σ, sim)
end

function Dynamics.set_quantum_derivative!(dσ, v, σ, sim::AbstractSimulation{<:AbstractEhrenfest})
    V = DynamicsUtils.calculate_density_propagator!(sim, v)
    DynamicsUtils.commutator!(dσ, V, σ, sim.calculator.tmp_mat_complex1)
    lmul!(-im, dσ)
end

include("ehrenfest.jl")
include("ehrenfest_rpmd.jl")

end # module
