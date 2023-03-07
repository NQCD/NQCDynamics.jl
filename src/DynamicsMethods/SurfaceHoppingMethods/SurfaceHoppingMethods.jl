
"""
    SurfaceHoppingMethods

Implementation for surface hopping methods.
"""
module SurfaceHoppingMethods

using DEDataArrays: DEDataArrays
using ComponentArrays: ComponentVector
using DiffEqBase: DiffEqBase
using LinearAlgebra: LinearAlgebra, lmul!
using OrdinaryDiffEq: OrdinaryDiffEq

using NQCDynamics:
    NQCDynamics,
    AbstractSimulation,
    Simulation,
    Calculators,
    DynamicsMethods,
    DynamicsUtils,
    Estimators,
    ndofs
using NQCModels: NQCModels, Model
using NQCBase: Atoms

mutable struct SurfaceHoppingVariables{T,A,Axes,S} <: DEDataArrays.DEDataVector{T}
    x::ComponentVector{T,A,Axes}
    state::S
end

DynamicsUtils.get_velocities(u::SurfaceHoppingVariables) = DynamicsUtils.get_velocities(u.x)
DynamicsUtils.get_positions(u::SurfaceHoppingVariables) = DynamicsUtils.get_positions(u.x)
DynamicsUtils.get_quantum_subsystem(u::SurfaceHoppingVariables) = DynamicsUtils.get_quantum_subsystem(u.x)

"""
Abstract type for all surface hopping methods.

Surface hopping methods follow the structure set out in this file.
The nuclear and electronic variables are propagated by the `motion!` function.
The surface hopping procedure is handled by the `HoppingCallback` which
uses the functions `check_hop!` and `execute_hop!` as its `condition` and `affect!`.

To add a new surface hopping scheme, you must create a new struct
and define methods for `evaluate_hopping_probability!`, `select_new_state`,
and `rescale_velocity!`.

See `fssh.jl` for an example implementation.
"""
abstract type SurfaceHopping <: DynamicsMethods.Method end

function DynamicsMethods.motion!(du, u, sim::Simulation{<:SurfaceHopping}, t)
    dr = DynamicsUtils.get_positions(du)
    dv = DynamicsUtils.get_velocities(du)
    dσ = DynamicsUtils.get_quantum_subsystem(du)

    r = DynamicsUtils.get_positions(u)
    v = DynamicsUtils.get_velocities(u)
    σ = DynamicsUtils.get_quantum_subsystem(u)

    set_state!(u, sim.method.state) # Make sure the state variables match, 

    DynamicsUtils.velocity!(dr, v, r, sim, t) # Set the velocity
    Calculators.update_electronics!(sim.calculator, r) # Calculate electronic quantities
    DynamicsUtils.acceleration!(dv, v, r, sim, t, sim.method.state) # Set the acceleration
    DynamicsUtils.set_quantum_derivative!(dσ, u, sim)
end

function DynamicsUtils.set_quantum_derivative!(dσ, u, sim::AbstractSimulation{<:SurfaceHopping})
    v = DynamicsUtils.get_hopping_velocity(sim, DynamicsUtils.get_velocities(u))
    σ = DynamicsUtils.get_quantum_subsystem(u)
    V = DynamicsUtils.calculate_density_matrix_propagator!(sim, v)
    DynamicsUtils.commutator!(dσ, V, σ)
    lmul!(-im, dσ)
end

function DynamicsMethods.create_problem(u0, tspan, sim::AbstractSimulation{<:SurfaceHopping})
    set_state!(sim.method, u0.state)
    OrdinaryDiffEq.ODEProblem(DynamicsMethods.motion!, u0, tspan, sim;
        callback=DynamicsMethods.get_callbacks(sim))
end

include("surface_hopping.jl")
include("fssh.jl")
include("iesh.jl")
include("rpsh.jl")
include("rpiesh.jl")
include("cme.jl")
export CME, BCME
include("rpcme.jl")

end # module
