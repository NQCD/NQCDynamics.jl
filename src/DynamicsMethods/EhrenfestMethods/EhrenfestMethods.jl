
module EhrenfestMethods

using LinearAlgebra: lmul!, eigvecs, diag, dot

using ....NonadiabaticMolecularDynamics:
    NonadiabaticMolecularDynamics,
    AbstractSimulation,
    Simulation,
    RingPolymerSimulation,
    Calculators,
    Dynamics
using ..Dynamics: DynamicsUtils

"""
Abstract type for Ehrenfest method.
"""
abstract type AbstractEhrenfest <: Dynamics.Method end

function Dynamics.motion!(du, u, sim::AbstractSimulation{<:AbstractEhrenfest}, t)
    dr = Dynamics.get_positions(du)
    dv = Dynamics.get_velocities(du)
    dσ = Dynamics.get_quantum_subsystem(du)

    r = Dynamics.get_positions(u)
    v = Dynamics.get_velocities(u)
    σ = Dynamics.get_quantum_subsystem(u)

    DynamicsUtils.velocity!(dr, v, r, sim, t)
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
