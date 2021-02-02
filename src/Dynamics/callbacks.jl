using DiffEqCallbacks

export CellBoundaryCallback

save_energy(u, t, integrator) =
    Models.energy(integrator.p.calculator.model, get_positions(u))

function create_energy_saving_callback()
    saved_values = SavedValues(Float64, Float64)
    SavingCallback(save_energy, saved_values), saved_values
end

outside_cell(u,t,integrator) = !check_atoms_in_cell(integrator.p.cell, get_positions(u))
function enforce_periodicity!(integrator)
    apply_cell_boundaries!(integrator.p.cell, get_positions(integrator.u))
end

const CellBoundaryCallback = DiscreteCallback(outside_cell, enforce_periodicity!)
