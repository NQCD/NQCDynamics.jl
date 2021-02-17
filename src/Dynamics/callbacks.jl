using DiffEqCallbacks

export CellBoundaryCallback
export create_energy_saving_callback
export create_terminating_callback

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

"""
    create_terminating_callback(func::Function)

Provide a function that returns true when the simulation should terminate.
"""
function create_terminating_callback(func::Function)
    DiscreteCallback(func, terminate!)
end
