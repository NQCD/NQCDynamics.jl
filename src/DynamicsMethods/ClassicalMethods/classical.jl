using OrdinaryDiffEq: OrdinaryDiffEq

"""
    Classical <: DynamicsMethods.Method

Type for performing classical molecular dynamics.

```jldoctest
sim = Simulation{Classical}(Atoms(:H), Harmonic())

# output

Simulation{Classical}:
  Atoms{1, Float64}([:H], UInt8[0x01], [1837.4715941070515])
  Harmonic{Float64, Float64, Float64}
  m: Float64 1.0
  ω: Float64 1.0
  r₀: Float64 0.0
  dofs: Int64 1
```
"""
struct Classical <: DynamicsMethods.Method end

function NQCDynamics.Simulation{Classical}(atoms::Atoms, model::Model; kwargs...)
    NQCDynamics.Simulation(atoms, model, Classical(); kwargs...)
end
function NQCDynamics.RingPolymerSimulation{Classical}(atoms::Atoms, model::Model, n_beads::Integer; kwargs...)
    NQCDynamics.RingPolymerSimulation(atoms, model, Classical(), n_beads::Integer; kwargs...)
end

function NQCDynamics.Simulation(atoms::Atoms, model::Model; kwargs...)
    NQCDynamics.Simulation{Classical}(atoms, model; kwargs...)
end
function NQCDynamics.RingPolymerSimulation(atoms::Atoms, model::Model, n_beads::Integer; kwargs...)
    NQCDynamics.RingPolymerSimulation{Classical}(atoms, model, n_beads; kwargs...)
end

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
