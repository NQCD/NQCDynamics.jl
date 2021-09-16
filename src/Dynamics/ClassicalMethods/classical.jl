
"""
A singleton type that simply labels the parent `AbstractSimulation` as classical.
"""
struct Classical <: Dynamics.Method end

"""
    motion!(du, u, sim::AbstractSimulation{<:Classical}, t)
    
Sets the time derivative for the positions and momenta contained within `u`.
"""
function Dynamics.motion!(du, u, sim::AbstractSimulation{<:Classical}, t)
    dr = get_positions(du)
    dv = get_velocities(du)
    r = get_positions(u)
    v = get_velocities(u)
    DynamicsUtils.velocity!(dr, v, r, sim, t)
    acceleration!(dv, v, r, sim, t)
end

function Dynamics.motion!(du, u, sim::RingPolymerSimulation{<:Classical}, t)
    dr = get_positions(du)
    dv = get_velocities(du)
    r = get_positions(u)
    v = get_velocities(u)
    DynamicsUtils.velocity!(dr, v, r, sim, t)
    ring_polymer_acceleration!(dv, v, r, sim, t)
end

"""
`f1` in `DifferentialEquations.jl` docs.
"""
function acceleration!(dv, v, r, sim::AbstractSimulation, t)
    Calculators.evaluate_derivative!(sim.calculator, r)
    dv .= -sim.calculator.derivative
    Dynamics.divide_by_mass!(dv, sim.atoms.masses)
end

function ring_polymer_acceleration!(dv, v, r, sim::RingPolymerSimulation, t)
    Calculators.evaluate_derivative!(sim.calculator, r)
    dv .= -sim.calculator.derivative
    Dynamics.divide_by_mass!(dv, sim.atoms.masses)
    Dynamics.apply_interbead_coupling!(dv, r, sim)
end

function Dynamics.create_problem(u0, tspan::Tuple, sim::Simulation{<:Classical})
    DynamicalODEProblem(acceleration!, velocity!, get_velocities(u0), get_positions(u0), tspan, sim)
end

function Dynamics.create_problem(u0, tspan::Tuple, sim::RingPolymerSimulation{<:Classical})
    DynamicalODEProblem(ring_polymer_acceleration!, velocity!, get_velocities(u0), get_positions(u0), tspan, sim)
end

Dynamics.select_algorithm(::AbstractSimulation{<:Classical}) = OrdinaryDiffEq.VelocityVerlet()
