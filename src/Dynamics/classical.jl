export Classical

"""
$(TYPEDEF)

A singleton type that simply labels the parent `AbstractSimulation` as classical.
"""
struct Classical <: Method end

"""
    motion!(du::DynamicalVariables, u::DynamicalVariables, sim::AbstractSimulation, t)
    
Sets the time derivative for the positions and momenta contained within `u`.

This is defined for the abstract types and acts as a fallback for all other dynamics methods.
"""
function motion!(du::DynamicalVariables, u::DynamicalVariables, sim::Simulation, t)
    dr = get_positions(du)
    dv = get_velocities(du)
    r = get_positions(u)
    v = get_velocities(u)
    velocity!(dr, v, r, sim, t)
    acceleration!(dv, v, r, sim, t)
end

function motion!(du::DynamicalVariables, u::DynamicalVariables, sim::RingPolymerSimulation, t)
    dr = get_positions(du)
    dv = get_velocities(du)
    r = get_positions(u)
    v = get_velocities(u)
    velocity!(dr, v, r, sim, t)
    acceleration!(dv, v, r, sim, t)
    apply_interbead_coupling!(dv, r, sim)
end

"""
`f2` in `DifferentialEquations.jl` docs.
"""
velocity!(dr, v, r, sim, t) = dr .= v

"""
`f1` in `DifferentialEquations.jl` docs.
"""
function acceleration!(dv, v, r, sim::AbstractSimulation, t)
    Calculators.evaluate_derivative!(sim.calculator, r)
    dv .= -sim.calculator.derivative ./ sim.atoms.masses'
end

"""
    apply_interbead_coupling!(du::DynamicalVariables, u::DynamicalVariables,
                              sim::RingPolymerSimulation)
    
Applies the force that arises from the harmonic springs between adjacent beads.

Only applies the force for atoms labelled as quantum within the `RingPolymerParameters`.
"""
function apply_interbead_coupling!(dr::RingPolymerArray, r::RingPolymerArray, sim::RingPolymerSimulation)
    for i in sim.beads.quantum_atoms
        for j=1:sim.DoFs
            dr[j,i,:] .-= 2sim.beads.springs*r[j,i,:]
        end
    end
end

function create_problem(u0::ClassicalDynamicals, tspan::Tuple, sim::AbstractSimulation{<:Classical})
    DynamicalODEProblem(acceleration!, velocity!, get_velocities(u0), get_positions(u0), tspan, sim)
end

function create_problem(u0::ArrayPartition, tspan::Tuple, sim::AbstractSimulation{<:Classical})
    DynamicalODEProblem(acceleration!, velocity!, u0.x[1], u0.x[2], tspan, sim)
end

select_algorithm(::AbstractSimulation{<:Classical}) = VelocityVerlet()
