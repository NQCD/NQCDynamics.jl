using DiffEqCallbacks
using Unitful
using UnitfulAtomic

export CellBoundaryCallback
export TerminatingCallback

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
