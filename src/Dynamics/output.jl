export create_saving_callback

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
potential(u, t, integrator) = Models.energy(integrator.p.calculator.model, get_positions(u))
energy(u, t, integrator) = evaluate_hamiltonian(integrator.p, u)
kinetic(u, t, integrator) = evaluate_kinetic_energy(integrator.p.atoms.masses, get_velocities(u))
u(u, t, integrator) = copy(u)
density_matrix(u, t, integrator) = copy(get_density_matrix(u))
state(u, t, integrator) = copy(u.state)
noise(u, t, integrator) = copy(integrator.W.dW) / sqrt(integrator.dt)
population(u, t, integrator) = get_population(integrator.p, u)
function friction(u, t, integrator)
    integrator.g(integrator.cache.gtmp,get_positions(u),integrator.p,t)
    copy(integrator.cache.gtmp)
end
