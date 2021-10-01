
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
    saved_values = SavedValues(Float64, NamedTuple{function_names})
    saving_function = get_saving_function(function_names)
    SavingCallback(saving_function, saved_values; saveat=saveat), saved_values
end
create_saving_callback(quantities::Symbol; saveat=[]) =
    create_saving_callback((quantities,), saveat=saveat)

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

"""
    get_saving_function(function_names::NTuple{N, Symbol})::Function where {N}

Take a Tuple of symbols representing the names of functions within this module
and return a function that will evaluate all of the provided functions.
"""
function get_saving_function(function_names::NTuple{N, Symbol})::Function where {N}

    output_functions = (getfield(DynamicsOutputs, f) for f in function_names)

    function saving(u, t, integrator)
        output_results = (f(u, t, integrator) for f in output_functions)
        NamedTuple{function_names}(output_results)
    end
end

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

"Evaluate the adiabatc population"
adiabatic_population(u, t, integrator) = Estimators.adiabatic_population(integrator.p, u)

"Evaluate the friction. This is used for MDEF only."
function friction(u, t, integrator)
    integrator.g(integrator.cache.gtmp,get_positions(u),integrator.p,t)
    copy(integrator.cache.gtmp)
end

end # module
