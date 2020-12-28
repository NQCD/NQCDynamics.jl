export Classical

"""
    Classical <: Method
    
A singleton type that simply labels the parent `AbstractSimulation` as classical.
"""
struct Classical <: Method end

"""
    motion!(du::DynamicalVariables, u::DynamicalVariables, sim::AbstractSimulation, t)
    
Sets the time derivative for the positions and momenta contained within `u`.

This is defined for the abstract types and acts as a fallback for all other dynamics methods.
"""
function motion!(du::DynamicalVariables, u::DynamicalVariables, sim::AbstractSimulation, t)
    set_velocity!(du, u, sim)
    set_force!(du, u, sim)
end

@doc raw"""
    set_velocity!(du::DynamicalVariables, u::DynamicalVariables, sim::AbstractSimulation)
    
Sets the velocity as ``\dot{R} = \frac{P}{M}``.
"""
function set_velocity!(du::DynamicalVariables, u::DynamicalVariables, sim::AbstractSimulation)
    get_positions(du) .= get_momenta(u) ./ sim.atoms.masses'
end

@doc raw"""
    set_force!(du::DynamicalVariables, u::DynamicalVariables, sim::AbstractSimulation) 
    
Sets the force as ``\dot{P_i} = -\frac{\partial V(R)}{\partial R_i}``
"""
function set_force!(du::DynamicalVariables, u::DynamicalVariables, sim::AbstractSimulation) 
    Calculators.evaluate_derivative!(sim.calculator, get_positions(u))
    get_momenta(du) .= -sim.calculator.derivative
end

"""
    motion!(du::DynamicalVariables, u::DynamicalVariables, sim::RingPolymerSimulation, T)
    
Sets the equations of motion for the ring polymer system.

Same as the classical equations but additionally includes the interbead coupling.
"""
function motion!(du::DynamicalVariables, u::DynamicalVariables, sim::RingPolymerSimulation, T)
    set_velocity!(du, u, sim)
    set_force!(du, u, sim)
    apply_interbead_coupling!(du, u, sim)
end

"""
    apply_interbead_coupling!(du::DynamicalVariables, u::DynamicalVariables,
                              sim::RingPolymerSimulation)
    
Applies the force that arises from the harmonic springs between adjacent beads.

Only applies the force for atoms labelled as quantum within the `RingPolymerParameters`.
"""
function apply_interbead_coupling!(du::DynamicalVariables, u::DynamicalVariables, sim::RingPolymerSimulation)
    for i in sim.beads.quantum_atoms
        for j=1:sim.DoFs
            get_momenta(du)[j,i,:] .-= 2sim.beads.springs*get_positions(u)[j,i,:] .* sim.atoms.masses[i]
        end
    end
end
