export Classical

"""
$(TYPEDEF)

A singleton type that simply labels the parent `AbstractSimulation` as classical.
"""
struct Classical <: Method end

"""
    motion!(du, u, sim::AbstractSimulation{<:Classical}, t)
    
Sets the time derivative for the positions and momenta contained within `u`.
"""
function motion!(du, u, sim::AbstractSimulation{<:Classical}, t)
    dr = get_positions(du)
    dv = get_velocities(du)
    r = get_positions(u)
    v = get_velocities(u)
    velocity!(dr, v, r, sim, t)
    acceleration!(dv, v, r, sim, t)
end

function motion!(du, u, sim::RingPolymerSimulation{<:Classical}, t)
    dr = get_positions(du)
    dv = get_velocities(du)
    r = get_positions(u)
    v = get_velocities(u)
    velocity!(dr, v, r, sim, t)
    ring_polymer_acceleration!(dv, v, r, sim, t)
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
    dv .= -sim.calculator.derivative
    divide_by_mass!(dv, sim.atoms.masses)
end

divide_by_mass!(dv, masses) = dv ./= masses'

function ring_polymer_acceleration!(dv, v, r, sim::RingPolymerSimulation, t)
    Calculators.evaluate_derivative!(sim.calculator, r)
    dv .= -sim.calculator.derivative
    divide_by_mass!(dv, sim.atoms.masses)
    apply_interbead_coupling!(dv, r, sim)
end

"""
    apply_interbead_coupling!(du::DynamicalVariables, u::DynamicalVariables,
                              sim::RingPolymerSimulation)
    
Applies the force that arises from the harmonic springs between adjacent beads.

Only applies the force for atoms labelled as quantum within the `RingPolymerParameters`.
"""
function apply_interbead_coupling!(dr::AbstractArray{T,3}, r::AbstractArray{T,3}, sim::RingPolymerSimulation) where {T}
    for i in axes(dr, 3)
        iplus = mod1(i+1, sim.beads.n_beads)
        iminus = mod1(i-1, sim.beads.n_beads)
        for j in sim.beads.quantum_atoms
            for k in axes(dr, 1)
                dr[k,j,i] -= sim.beads.ω_n² * (2r[k,j,i] - r[k,j,iplus] - r[k,j,iminus])
            end
        end
    end
    return nothing
end

function create_problem(u0, tspan::Tuple, sim::Simulation{<:Classical})
    DynamicalODEProblem(acceleration!, velocity!, get_velocities(u0), get_positions(u0), tspan, sim)
end

function create_problem(u0, tspan::Tuple, sim::RingPolymerSimulation{<:Classical})
    DynamicalODEProblem(ring_polymer_acceleration!, velocity!, get_velocities(u0), get_positions(u0), tspan, sim)
end

select_algorithm(::AbstractSimulation{<:Classical}) = VelocityVerlet()
