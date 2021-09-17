using RecursiveArrayTools: ArrayPartition
using DiffEqCallbacks: SavedValues, SavingCallback

using NonadiabaticModels: NonadiabaticModels
using NonadiabaticMolecularDynamics:
    evaluate_hamiltonian,
    evaluate_kinetic_energy,
    get_centroid

using ..DynamicsUtils:
    get_positions,
    get_velocities,
    get_quantum_subsystem

using NonadiabaticMolecularDynamics.Estimators:
    diabatic_population,
    adiabatic_population

create_saving_callback(quantities::Symbol; saveat=[]) = create_saving_callback((quantities,), saveat=saveat)

function create_saving_callback(quantities::NTuple{N, Symbol}; saveat=[]) where {N}
    saved_values = SavedValues(Float64, NamedTuple{quantities})
    saving_function = get_saving_function(NamedTuple{quantities})
    SavingCallback(saving_function, saved_values; saveat=saveat), saved_values
end

function get_saving_function(::Type{savevalType})::Function where {savevalType}

    evaluate_field(field, u, t, integrator) = @eval $field($u, $t, $integrator)

    function saving(u, t, integrator)::savevalType
        output = [evaluate_field(field, u, t, integrator) for field in fieldnames(savevalType)]
        savevalType(output)
    end
end

force(u, t, integrator) = -copy(integrator.p.calculator.derivative)
velocity(u, t, integrator) = copy(get_velocities(u))
position(u, t, integrator) = copy(get_positions(u))
centroid_position(u, t, integrator) = get_centroid(get_positions(u))
potential(u, t, integrator) = NonadiabaticModels.potential(integrator.p.calculator.model, get_positions(u))[1]
hamiltonian(u, t, integrator) = evaluate_hamiltonian(integrator.p, u)
kinetic(u, t, integrator) = evaluate_kinetic_energy(integrator.p.atoms.masses, get_velocities(u))
u(u, t, integrator) = copy(u)
u(u::ArrayPartition, t, integrator) = ComponentVector(v=copy(get_velocities(u)), r=copy(get_positions(u)))
quantum_subsystem(u, t, integrator) = copy(get_quantum_subsystem(u))
state(u, t, integrator) = copy(u.state)
noise(u, t, integrator) = copy(integrator.W.dW) / sqrt(integrator.dt)
population(u, t, integrator) = diabatic_population(integrator.p, u)
adiabatic_population(u, t, integrator) = adiabatic_population(integrator.p, u)
function friction(u, t, integrator)
    integrator.g(integrator.cache.gtmp,get_positions(u),integrator.p,t)
    copy(integrator.cache.gtmp)
end
