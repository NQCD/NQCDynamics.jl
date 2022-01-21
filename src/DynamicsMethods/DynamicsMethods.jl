"""
This module contains functions and types necessary for performing
nonadiabatic molecular dynamics.

Dynamics is performed using [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/).
As such, this module is centered around the implementation of the functions
necessary to integrate the dynamics.

For deterministic Hamiltonian methods, the central function is [`DynamicsMethods.motion!`](@ref),
which is the inplace form of the function
to be integrated by [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/).

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
using ComponentArrays: ComponentVector

using NonadiabaticMolecularDynamics: AbstractSimulation, Simulation, RingPolymerSimulation,
    DynamicsOutputs
using NonadiabaticDynamicsBase: austrip_kwargs

export DynamicsVariables
export run_trajectory

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

"Provides the DEProblem for each type of simulation."
create_problem(u0, tspan, sim) =
    OrdinaryDiffEq.ODEProblem(motion!, u0, tspan, sim; callback=get_callbacks(sim))

"Choose a default algorithm for solving the differential equation."
select_algorithm(::AbstractSimulation) = OrdinaryDiffEq.VCABM5()

"Select the default callbacks for this simulation type."
get_callbacks(::AbstractSimulation) = nothing

"""
    DynamicsVariables(::AbstractSimulation, args...)

For each dynamics method this function is implemented to provide the variables for the
dynamics in the appropriate format.

By default, `DynamicsVariables` is set up for the classical case and takes
`sim`, `v`, `r` as arguments and returns a `ComponentVector(v=v, r=r)`
which is used as a container for the velocities and positions during
classical dynamics.
"""
DynamicsVariables(::AbstractSimulation, v, r) = ComponentVector(v=v, r=r)

"""
    run_trajectory(u0, tspan::Tuple, sim::AbstractSimulation;
                   output=(:u,), saveat=[], algorithm=select_algorithm(sim), kwargs...)

Solve a single trajectory starting from `u0` with a timespan `tspan`
for the simulation `sim`.

# Keyword arguments

`output` specifies the quantities that should be saved during the dynamics simulation.
The options for this keyword are any of the functions found in
`src/DynamicsMethods/output.jl`.

The rest of the keywords are the usual arguments for the `solve` function from
`DifferentialEquations.jl`.
It is possible to use `Unitful` quantities for any of the arguments since these are
automatically converted to atomic units internally.

# Output

The function returns a `Table` from `TypedTables.jl` with columns for time `t` and
every quantity specified in the `output` tuple.
"""
function run_trajectory(u0, tspan::Tuple, sim::AbstractSimulation;
        output=(:u,), saveat=[], callback=nothing, algorithm=select_algorithm(sim),
        kwargs...)

    if !(output isa Tuple)
        output = (output,)
    end

    problem = create_problem(u0, austrip.(tspan), sim)

    saving_callback, vals = DynamicsOutputs.create_saving_callback(output; saveat=austrip.(saveat))
    callback_set = DiffEqBase.CallbackSet(callback, saving_callback, get_callbacks(sim))
    problem = DiffEqBase.remake(problem, callback=callback_set)

    stripped_kwargs = austrip_kwargs(;kwargs...)
    DiffEqBase.solve(problem, algorithm; stripped_kwargs...)
    out = [(;zip(output, val)...) for val in vals.saveval]
    TypedTables.Table(t=vals.t, out)
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