using DiffEqBase: DiscreteCallback

save_impurity(u, t, integrator) =
    Models.impurity_summary(integrator.p.calculator.model, get_positions(u), u.state, u.x.x[3])

outside_cell(u,t,integrator) = !check_atoms_in_cell(integrator.p.cell, get_positions(u))
function enforce_periodicity!(integrator)
    apply_cell_boundaries!(integrator.p.cell, get_positions(integrator.u))
end

const CellBoundaryCallback = DiscreteCallback(outside_cell, enforce_periodicity!)

"""
    TerminatingCallback(func::Function)

Provide a function that returns true when the simulation should terminate.
"""
TerminatingCallback(func::Function) = DiscreteCallback(func, terminate!)
