using OrdinaryDiffEq: OrdinaryDiffEq

"""
A singleton type that simply labels the parent `AbstractSimulation` as classical.
"""
struct Classical <: Dynamics.Method end

"""
    motion!(du, u, sim::AbstractSimulation{<:Classical}, t)
    
Sets the time derivative for the positions and momenta contained within `u`.
"""
function Dynamics.motion!(du, u, sim::AbstractSimulation{<:Classical}, t)
    dr = Dynamics.get_positions(du)
    dv = Dynamics.get_velocities(du)
    r = Dynamics.get_positions(u)
    v = Dynamics.get_velocities(u)
    DynamicsUtils.velocity!(dr, v, r, sim, t)
    acceleration!(dv, v, r, sim, t)
end

function Dynamics.motion!(du, u, sim::RingPolymerSimulation{<:Classical}, t)
    dr = Dynamics.get_positions(du)
    dv = Dynamics.get_velocities(du)
    r = Dynamics.get_positions(u)
    v = Dynamics.get_velocities(u)
    DynamicsUtils.velocity!(dr, v, r, sim, t)
    ring_polymer_acceleration!(dv, v, r, sim, t)
end

"""
`f1` in `DifferentialEquations.jl` docs.
"""
function acceleration!(dv, v, r, sim::AbstractSimulation, t)
    Calculators.evaluate_derivative!(sim.calculator, r)
    dv .= -sim.calculator.derivative
    DynamicsUtils.divide_by_mass!(dv, sim.atoms.masses)
end

function ring_polymer_acceleration!(dv, v, r, sim::RingPolymerSimulation, t)
    Calculators.evaluate_derivative!(sim.calculator, r)
    dv .= -sim.calculator.derivative
    DynamicsUtils.divide_by_mass!(dv, sim.atoms.masses)
    DynamicsUtils.apply_interbead_coupling!(dv, r, sim)
end

function Dynamics.create_problem(u0, tspan::Tuple, sim::Simulation{<:Classical})
    OrdinaryDiffEq.DynamicalODEProblem(acceleration!, DynamicsUtils.velocity!,
        Dynamics.get_velocities(u0), Dynamics.get_positions(u0), tspan, sim)
end

function Dynamics.create_problem(u0, tspan::Tuple, sim::RingPolymerSimulation{<:Classical})
    OrdinaryDiffEq.DynamicalODEProblem(ring_polymer_acceleration!, DynamicsUtils.velocity!,
        Dynamics.get_velocities(u0), Dynamics.get_positions(u0), tspan, sim)
end

Dynamics.select_algorithm(::AbstractSimulation{<:Classical}) = OrdinaryDiffEq.VelocityVerlet()
