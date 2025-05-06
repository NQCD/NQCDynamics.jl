
module EhrenfestMethods

using LinearAlgebra: lmul!, eigvecs, diag, dot

using NQCDynamics:
    NQCDynamics,
    AbstractSimulation,
    Simulation,
    RingPolymerSimulation,
    DynamicsMethods,
    DynamicsUtils,
    Estimators
using NQCCalculators
using NQCModels: NQCModels, Model
using NQCBase: Atoms

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
    NQCCalculators.update_electronics!(sim.cache, r)
    DynamicsUtils.acceleration!(dv, v, r, sim, t, σ)
    DynamicsUtils.set_quantum_derivative!(dσ, u, sim)
end

function DynamicsUtils.set_quantum_derivative!(dσ, u, sim::AbstractSimulation{<:AbstractEhrenfest})
    v = DynamicsUtils.get_hopping_velocity(sim, DynamicsUtils.get_velocities(u))
    σ = DynamicsUtils.get_quantum_subsystem(u)
    r = DynamicsUtils.get_positions(u)
    eigenvalues = DynamicsUtils.get_hopping_eigenvalues(sim, r)
    propagator = sim.method.density_propagator
    d = DynamicsUtils.get_hopping_nonadiabatic_coupling(sim, r)
    V = DynamicsUtils.calculate_density_matrix_propagator!(propagator, v, d, eigenvalues)

    tmp1 = sim.method.tmp_complex_matrix
    tmp2 = sim.method.tmp_complex_matrix2
    copy!(tmp2, σ) # Copy to complex matrix for faster mul!
    DynamicsUtils.commutator!(tmp1, V, tmp2)
    copy!(dσ, tmp1)
    lmul!(-im, dσ)
end

include("ehrenfest.jl")
include("ehrenfest_rpmd.jl")
include("ehrenfest_na.jl")
export EhrenfestNA
include("rpehrenfest_na.jl")

end # module
