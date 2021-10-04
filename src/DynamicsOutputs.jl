
module DynamicsOutputs

using RecursiveArrayTools: ArrayPartition
using DiffEqCallbacks: SavedValues, SavingCallback
using ComponentArrays: ComponentVector

using NonadiabaticModels: NonadiabaticModels
using NonadiabaticMolecularDynamics:
    Estimators,
    DynamicsUtils

using ..DynamicsUtils:
    get_positions,
    get_velocities,
    get_quantum_subsystem

using ..RingPolymers: get_centroid

"""
    create_saving_callback(quantities::NTuple{N, Symbol}; saveat=[]) where {N}

Get the `SavingCallback` that will populate `saved_values` with the result obtained
by evaluating the functions provided in `function_names`.
"""
function create_saving_callback(function_names::NTuple{N, Symbol}; saveat=[]) where {N}
    saved_values = SavedValues(Float64, Vector{Any})
    saving_function = OutputSaver(function_names)
    SavingCallback(saving_function, saved_values; saveat=saveat), saved_values
end

struct OutputSaver{N,F}
    function_names::NTuple{N, Symbol}
    output_functions::F
end

"""
    OutputSaver(function_names::NTuple{N, Symbol}) where {N}

Used to obtain the outputs for all functions given in `function_names`.
"""
function OutputSaver(function_names::NTuple{N, Symbol}) where {N}
    output_functions = tuple([getfield(DynamicsOutputs, f) for f in function_names]...)
    OutputSaver(function_names, output_functions)
end

"""
    (output::OutputSaver)(u, t, integrator)

Evaluates every function listed in `output.function_names` and returns all the results.
"""
function (output::OutputSaver)(u, t, integrator)
    out = Any[]
    sizehint!(out, length(output.function_names))
    evaluate_output!(out, u, t, integrator, output.output_functions...)
    return out
end

"""
Used to recursively evaluate every function for the [`OutputSaver`](@ref). 

See here for a description of why it is written like this: 
https://stackoverflow.com/questions/55840333/type-stability-for-lists-of-closures
"""
function evaluate_output!(out, u, t, integrator, f::F, output_functions...) where {F}
    push!(out, f(u, t, integrator))
    return evaluate_output!(out, u, t, integrator, output_functions...)
end
evaluate_output!(out, u, t, integrator) = out

export force
export velocity
export position
export centroid_position
export potential
export hamiltonian
export kinetic
export u
export quantum_subsystem
export state
export noise
export population
export adiabatic_population
export friction

"Get the force"
force(u, t, integrator) = -copy(integrator.p.calculator.derivative)

"Get the velocity"
velocity(u, t, integrator) = copy(get_velocities(u))

"Get the position"
position(u, t, integrator) = copy(get_positions(u))

"Get the position of the ring polymer centroid"
centroid_position(u, t, integrator) = get_centroid(get_positions(u))

"Evaluate the potential from the model"
potential(u, t, integrator) = NonadiabaticModels.potential(integrator.p.calculator.model, get_positions(u))[1]

"Evaluate the classical Hamiltonian"
hamiltonian(u, t, integrator) = DynamicsUtils.classical_hamiltonian(integrator.p, u)

"Evaluate the classical kinetic energy"
kinetic(u, t, integrator) = DynamicsUtils.classical_kinetic_energy(integrator.p, get_velocities(u))

"Get all the dynamics variables. This is the default"
u(u, t, integrator) = copy(u)
u(u::ArrayPartition, t, integrator) = ComponentVector(v=copy(get_velocities(u)), r=copy(get_positions(u)))

"""
Get the quantum subsystem of the dynamics variables.
Requires that `DynamicsUtils.get_quantum_subsystem` is implemented for your chosen method.
"""
quantum_subsystem(u, t, integrator) = copy(get_quantum_subsystem(u))

"""
Get the currently occupied state from the dynamics variables.
Requires that the dynamics variable has a field `state`.
Currently this is for surface hopping methods only.
Use [`population`](@ref) or [`adiabatic_population`](@ref) for most other methods.
"""
state(u, t, integrator) = copy(u.state)

"Get the noise along the stochastic trajectory"

noise(u, t, integrator) = copy(integrator.W.dW) / sqrt(integrator.dt)

"Evaluate the diabatic population"
population(u, t, integrator) = Estimators.diabatic_population(integrator.p, u)

"Evaluate the adiabatic population"
adiabatic_population(u, t, integrator) = Estimators.adiabatic_population(integrator.p, u)

"Evaluate the friction. This is used for MDEF only."
function friction(u, t, integrator)
    integrator.g(integrator.cache.gtmp,get_positions(u),integrator.p,t)
    copy(integrator.cache.gtmp)
end

end # module
