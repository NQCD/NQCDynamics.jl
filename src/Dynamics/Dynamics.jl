"""
This module contains functions and types necessary for performing
nonadiabatic molecular dynamics.

Dynamics is performed using [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/).
As such, this module is centered around the implementation of the functions
necessary to integrate the dynamics.

For deterministic Hamiltonian methods, the central function is [`Dynamics.motion!`](@ref),
which is the inplace form of the function
to be integrated by [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/).

For stochastic methods, `motion!` provides the deterministic part of the equation,
and `random_force!` should be implemented to provide the stochastic part.

Further, methods that have discontinuities, such as surface hopping, use the
[callback interface](https://diffeq.sciml.ai/stable/features/callback_functions/#callbacks)
provided by `DifferentialEquations.jl`.
"""
module Dynamics

export motion!
export random_force!

using ..NonadiabaticMolecularDynamics
using DiffEqBase
using StochasticDiffEq
using OrdinaryDiffEq
using RecursiveArrayTools: ArrayPartition

"""
Each type of dynamics subtypes `Method` which is passed to
the `AbstractSimulation` as a parameter to determine the type of
dynamics desired.
"""
abstract type Method end

"""
    motion!(du, u, sim, t)
    
As per `DifferentialEquations.jl`, this function is implemented for
each method and defines the time derivatives of the `DynamicalVariables`.

We require that each implementation ensures `du` and `u` are subtypes
of `DynamicalVariables` and `sim` subtypes `AbstractSimulation`.
"""
function motion! end

"""
    random_force!(du, u, sim, t)
    
Similarly to [`Dynamics.motion!`](@ref), this function is directly passed
to an `SDEProblem` to integrate stochastic dynamics.
"""
function random_force! end

function run_trajectory(u0::DynamicalVariables, tspan::Tuple, sim::AbstractSimulation; kwargs...)
    solve(create_problem(u0, tspan, sim), select_algorithm(sim); kwargs...)
end

function create_problem(u0::DynamicalVariables, tspan::Tuple, sim::AbstractSimulation)
    ODEProblem(motion!, u0, tspan, sim)
end

select_algorithm(::AbstractSimulation) = Tsit5()

include("classical.jl")
include("langevin.jl")
include("mdef.jl")
include("fssh.jl")
include("iesh.jl")
include("fermionic_ring_polymer.jl")
include("nrpmd.jl")

include("ensembles.jl")

end # module