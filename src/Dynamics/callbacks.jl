using DiffEqCallbacks

save_energy(u, t, integrator) =
    Models.energy(integrator.p.calculator.model, get_positions(u))

function create_energy_saving_callback()
    saved_values = SavedValues(Float64, Float64)
    a = SavingCallback(save_energy, saved_values)
    SavingCallback(save_energy, saved_values), saved_values
end
