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
module DynamicsMethods

using DiffEqBase: DiffEqBase
using TypedTables: TypedTables
using UnitfulAtomic: austrip
using Reexport: @reexport
using OrdinaryDiffEq: OrdinaryDiffEq

using ..NonadiabaticMolecularDynamics: AbstractSimulation, Simulation, RingPolymerSimulation
using NonadiabaticDynamicsBase: austrip_kwargs

export DynamicsVariables

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
Provides the DEProblem for each type of simulation.
"""
create_problem(u0, tspan, sim) =
    ODEProblem(motion!, u0, tspan, sim; callback=get_callbacks(sim))

select_algorithm(::AbstractSimulation) = OrdinaryDiffEq.VCABM5()
get_callbacks(::AbstractSimulation) = nothing

function set_quantum_derivative! end

DynamicsVariables(::AbstractSimulation, v, r) = ComponentVector(v=v, r=r)

"""
    run_trajectory(u0::DynamicalVariables, tspan::Tuple, sim::AbstractSimulation;
        output=(:u,), callback=nothing, kwargs...)

Solve a single trajectory.
"""
function run_trajectory(u0, tspan::Tuple, sim::AbstractSimulation;
        output=(:u,), saveat=[], callback=nothing, algorithm=select_algorithm(sim),
        kwargs...)

    saving_callback, vals = DynamicsUtils.create_saving_callback(output; saveat=austrip.(saveat))
    callback_set = DiffEqBase.CallbackSet(callback, saving_callback, get_callbacks(sim))
    problem = create_problem(u0, austrip.(tspan), sim)
    problem = DiffEqBase.remake(problem, callback=callback_set)

    stripped_kwargs = austrip_kwargs(;kwargs...)
    DiffEqBase.solve(problem, algorithm; stripped_kwargs...)
    TypedTables.Table(t=vals.t, vals.saveval)
end

include("ClassicalMethods/ClassicalMethods.jl")
@reexport using .ClassicalMethods

include("SurfaceHoppingMethods/SurfaceHoppingMethods.jl")
@reexport using .SurfaceHoppingMethods

include("MappingVariableMethods/MappingVariableMethods.jl")
@reexport using .MappingVariableMethods

include("EhrenfestMethods/EhrenfestMethods.jl")
@reexport using .EhrenfestMethods

include("IntegrationAlgorithms/IntegrationAlgorithms.jl")

end # module