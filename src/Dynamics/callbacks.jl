using DiffEqCallbacks
using Unitful
using UnitfulAtomic

export CellBoundaryCallback
export create_terminating_callback
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

force(u, t, integrator) = -integrator.p.calculator.derivative
velocity(u, t, integrator) = copy(get_velocities(u))
position(u, t, integrator) = copy(get_positions(u))
potential(u, t, integrator) = Models.energy(integrator.p.calculator.model, get_positions(u))
energy(u, t, integrator) = evaluate_hamiltonian(integrator.p, u)
u(u, t, integrator) = copy(u)
density_matrix(u, t, integrator) = copy(get_density_matrix(u))
state(u, t, integrator) = copy(u.state)
noise(u, t, integrator) = copy(integrator.W.dW) / sqrt(integrator.dt)
function friction(u, t, integrator)
    integrator.g(integrator.cache.gtmp,get_positions(u),integrator.p,t)
    auconvert.(u"ps^-1", integrator.cache.gtmp)
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
