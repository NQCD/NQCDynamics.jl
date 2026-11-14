using DiffEqBase: DiscreteCallback, terminate!
using NQCBase: apply_cell_boundaries!, check_atoms_in_cell

outside_cell(u,t,integrator) = !check_atoms_in_cell(integrator.p.cell, get_positions(u))
function enforce_periodicity!(integrator)
    apply_cell_boundaries!(integrator.p.cell, get_positions(integrator.u))
end

"""
    CellBoundaryCallback()

Whenever atoms leave the simulation cell, enforce the periodicity by wrapping the positions
at the cell boundary.
"""
CellBoundaryCallback() = DiscreteCallback(outside_cell, enforce_periodicity!)

"""
    TerminateCellCallback()

If the atoms leave the simulation cell, terminate the simulation.
"""
TerminateCellCallback() = DiscreteCallback(outside_cell, terminate!)

"""
    TerminatingCallback(func)

Provide a function that returns true when the simulation should terminate.
"""
TerminatingCallback(func) = DiscreteCallback(func, terminate!)
