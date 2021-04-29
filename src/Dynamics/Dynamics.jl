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
using UnitfulAtomic, Unitful
using DocStringExtensions
using TypedTables

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

"""
    run_trajectory(u0::DynamicalVariables, tspan::Tuple, sim::AbstractSimulation;
        output=(:u,), callback=nothing, kwargs...)

Solve a single trajectory.
"""
function run_trajectory(u0::DynamicalVariables, tspan::Tuple, sim::AbstractSimulation; output=(:u,), saveat=[], callback=nothing, kwargs...)
    stripped_kwargs = austrip_kwargs(;kwargs...)
    saving_callback, vals = create_saving_callback(output; saveat=austrip.(saveat))
    callback_set = CallbackSet(callback, saving_callback, get_callbacks(sim))
    problem = create_problem(u0, austrip.(tspan), sim)
    problem = remake(problem, callback=callback_set)
    solve(problem, select_algorithm(sim); stripped_kwargs...)
    Table(t=auconvert.(u"fs", vals.t), vals.saveval)
end

"""
Provides the DEProblem for each type of simulation.
"""
create_problem(u0, tspan, sim) = ODEProblem(motion!, u0, tspan, sim)

select_algorithm(::AbstractSimulation) = Tsit5()
get_callbacks(::AbstractSimulation) = nothing

include("classical.jl")
include("langevin.jl")
include("mdef.jl")
include("SurfaceHopping/SurfaceHopping.jl")
include("fermionic_ring_polymer.jl")
include("nrpmd.jl")

include("algorithms/mdef_baoab.jl")
include("algorithms/bcocb.jl")

include("callbacks.jl")
include("output.jl")
include("plot.jl")

end # module