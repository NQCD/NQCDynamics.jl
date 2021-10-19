
module EhrenfestMethods

using LinearAlgebra: lmul!, eigvecs, diag, dot

using NonadiabaticMolecularDynamics:
    NonadiabaticMolecularDynamics,
    AbstractSimulation,
    Simulation,
    RingPolymerSimulation,
    Calculators,
    DynamicsMethods,
    DynamicsUtils,
    Estimators
using NonadiabaticModels: NonadiabaticModels, Model
using NonadiabaticDynamicsBase: Atoms

"""
Abstract type for Ehrenfest method.
"""
abstract type AbstractEhrenfest <: DynamicsMethods.Method end

function DynamicsMethods.motion!(du, u, sim::Simulation{<:AbstractEhrenfest}, t)
    dr = DynamicsUtils.get_positions(du)
    dv = DynamicsUtils.get_velocities(du)
    dσ = DynamicsUtils.get_quantum_subsystem(du)

    r = DynamicsUtils.get_positions(u)
    v = DynamicsUtils.get_velocities(u)
    σ = DynamicsUtils.get_quantum_subsystem(u)

    DynamicsUtils.velocity!(dr, v, r, sim, t)
    Calculators.update_electronics!(sim.calculator, r)
    acceleration!(dv, v, r, sim, t, σ)
    DynamicsUtils.set_quantum_derivative!(dσ, v, σ, sim)
end

function DynamicsUtils.set_quantum_derivative!(dσ, v, σ, sim::AbstractSimulation{<:AbstractEhrenfest})
    V = DynamicsUtils.calculate_density_matrix_propagator!(sim, v)
    DynamicsUtils.commutator!(dσ, V, σ, sim.calculator.tmp_mat_complex1)
    lmul!(-im, dσ)
end

include("ehrenfest.jl")
include("ehrenfest_rpmd.jl")

end # module
