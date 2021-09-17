using OrdinaryDiffEq: OrdinaryDiffEq

"""
A singleton type that simply labels the parent `AbstractSimulation` as classical.
"""
struct Classical <: DynamicsMethods.Method end

"""
    motion!(du, u, sim::AbstractSimulation{<:Classical}, t)
    
Sets the time derivative for the positions and momenta contained within `u`.
"""
function DynamicsMethods.motion!(du, u, sim::AbstractSimulation{<:Classical}, t)
    dr = DynamicsUtils.get_positions(du)
    dv = DynamicsUtils.get_velocities(du)
    r = DynamicsUtils.get_positions(u)
    v = DynamicsUtils.get_velocities(u)
    DynamicsUtils.velocity!(dr, v, r, sim, t)
    acceleration!(dv, v, r, sim, t)
end

function DynamicsMethods.motion!(du, u, sim::RingPolymerSimulation{<:Classical}, t)
    dr = DynamicsUtils.get_positions(du)
    dv = DynamicsUtils.get_velocities(du)
    r = DynamicsUtils.get_positions(u)
    v = DynamicsUtils.get_velocities(u)
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

function DynamicsMethods.create_problem(u0, tspan::Tuple, sim::Simulation{<:Classical})
    OrdinaryDiffEq.DynamicalODEProblem(acceleration!, DynamicsUtils.velocity!,
        DynamicsUtils.get_velocities(u0), DynamicsUtils.get_positions(u0), tspan, sim)
end

function DynamicsMethods.create_problem(u0, tspan::Tuple, sim::RingPolymerSimulation{<:Classical})
    OrdinaryDiffEq.DynamicalODEProblem(ring_polymer_acceleration!, DynamicsUtils.velocity!,
        DynamicsUtils.get_velocities(u0), DynamicsUtils.get_positions(u0), tspan, sim)
end

DynamicsMethods.select_algorithm(::AbstractSimulation{<:Classical}) = OrdinaryDiffEq.VelocityVerlet()
