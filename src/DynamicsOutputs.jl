
"""
    DynamicsOutputs

Infrastructure for saving quantities during trajectories.

Defines all available options that can be used within the `output` tuple along with
the functions that perform the saving operation.
"""
module DynamicsOutputs

using RecursiveArrayTools: ArrayPartition
using DiffEqCallbacks: SavedValues, SavingCallback
using ComponentArrays: ComponentVector

using NQCModels: NQCModels
using NQCDynamics:
    Estimators,
    DynamicsUtils,
    AbstractSimulation

using ..DynamicsUtils:
    get_positions,
    get_velocities,
    get_quantum_subsystem

using ..RingPolymers: get_centroid

struct EnsembleSaver{N,F,S<:AbstractSimulation}
    function_names::NTuple{N, Symbol}
    output_functions::F
    sim::S
end

function EnsembleSaver(function_names::NTuple{N, Symbol}, sim) where {N}
    output_functions = tuple([getfield(DynamicsOutputs, f) for f in function_names]...)
    EnsembleSaver(function_names, output_functions, sim)
end

function output_template(output::EnsembleSaver, savepoints, u0)
    sol = (t=savepoints, u=[u0 for _ in savepoints])
    out = output(sol, 0)[1]
    for i in eachindex(out)
        for j in eachindex(out[i])
            out[i][j] = zero(out[i][j])
        end
    end
    return out
end

function (output::EnsembleSaver)(sol, _)
    out = [Any[] for _ in 1:length(sol.t)]
    for (i, (t, u)) in enumerate(zip(sol.t, sol.u))
        sizehint!(out[i], length(output.function_names)+1)
        push!(out[i], t)
        evaluate_output!(out[i], u, t, output.sim, output.output_functions...)
    end
    return (out, false)
end

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
    evaluate_output!(out, u, t, integrator.p, output.output_functions...)
    return out
end

"""
Used to recursively evaluate every function for the [`OutputSaver`](@ref). 

See here for a description of why it is written like this: 
https://stackoverflow.com/questions/55840333/type-stability-for-lists-of-closures
"""
function evaluate_output!(out, u, t, sim, f::F, output_functions...) where {F}
    push!(out, f(u, t, sim))
    return evaluate_output!(out, u, t, sim, output_functions...)
end
evaluate_output!(out, u, t, sim) = out

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
force(u, t, sim) = -copy(sim.calculator.derivative)

"Get the velocity"
velocity(u, t, sim) = copy(get_velocities(u))

"Get the velocity of the ring polymer centroid"
centroid_velocity(u, t, sim) = get_centroid(get_velocities(u))

"Get the position"
position(u, t, sim) = copy(get_positions(u))

"Get the position of the ring polymer centroid"
centroid_position(u, t, sim) = get_centroid(get_positions(u))

"Evaluate the potential from the model"
potential(u, t, sim) = DynamicsUtils.classical_potential_energy(sim, u)

"Evaluate the classical Hamiltonian"
hamiltonian(u, t, sim) = DynamicsUtils.classical_hamiltonian(sim, u)

"Evaluate the classical kinetic energy"
kinetic(u, t, sim) = DynamicsUtils.classical_kinetic_energy(sim, get_velocities(u))

"Get all the dynamics variables. This is the default"
u(u, t, sim) = copy(u)
u(u::ArrayPartition, t, sim) = ComponentVector(v=copy(get_velocities(u)), r=copy(get_positions(u)))

"""
Get the quantum subsystem of the dynamics variables.
Requires that `DynamicsUtils.get_quantum_subsystem` is implemented for your chosen method.
"""
quantum_subsystem(u, t, sim) = copy(get_quantum_subsystem(u))

mapping_position(u, t, sim) = copy(DynamicsUtils.get_mapping_positions(u))
mapping_momentum(u, t, sim) = copy(DynamicsUtils.get_mapping_momenta(u))

"""
Get the currently occupied state from the dynamics variables.
Requires that the dynamics variable has a field `state`.
Currently this is for surface hopping methods only.
Use [`population`](@ref) or [`adiabatic_population`](@ref) for most other methods.
"""
state(u, t, sim) = copy(u.state)

"Evaluate the diabatic population"
population(u, t, sim) = Estimators.diabatic_population(sim, u)
total_population(u, t, sim) = sum(Estimators.diabatic_population(sim, u))

"Evaluate the adiabatic population"
adiabatic_population(u, t, sim) = Estimators.adiabatic_population(sim, u)

end # module
