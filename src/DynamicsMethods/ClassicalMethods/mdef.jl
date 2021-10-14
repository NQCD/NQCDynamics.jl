
using NonadiabaticMolecularDynamics.Calculators: DiabaticFrictionCalculator

abstract type AbstractMDEF <: DynamicsMethods.Method end

"""
```math
dr = v dt\\\\
dv = -\\Delta U/M dt - \\Gamma v dt + \\sigma \\sqrt{2\\Gamma} dW
```
``\\Gamma`` is the friction tensor with units of inverse time.
For thermal dynamics we set ``\\sigma = \\sqrt{kT / M}``,
where ``T`` is the electronic temperature.

This is integrated using the BAOAB algorithm where the friction "O" step is performed
in the tensor's eigenbasis. See `src/dynamics/mdef_baoab.jl` for details.
"""
struct MDEF{M} <: AbstractMDEF
    mass_scaling::M
end

MDEF(masses::AbstractVector, DoFs::Integer) = MDEF(get_mass_scale_matrix(masses, DoFs))

function get_mass_scale_matrix(masses::AbstractVector, DoFs::Integer)
    vect = repeat(sqrt.(masses), inner=DoFs)
    vect * vect'
end

function NonadiabaticMolecularDynamics.Simulation{MDEF}(atoms::Atoms, model::Model; kwargs...)
    NonadiabaticMolecularDynamics.Simulation(atoms, model, MDEF(atoms.masses, ndofs(model));  kwargs...)
end
function NonadiabaticMolecularDynamics.RingPolymerSimulation{MDEF}(atoms::Atoms, model::Model, n_beads::Integer; kwargs...)
    NonadiabaticMolecularDynamics.RingPolymerSimulation(atoms, model, MDEF(atoms.masses, ndofs(model)), n_beads; kwargs...)
end

"""
    acceleration!(dv, v, r, sim::Simulation{MDEF,<:DiabaticFrictionCalculator}, t)

Sets acceleration due to ground state force when using a `DiabaticFrictionModel`.
"""
function acceleration!(dv, v, r, sim::Simulation{MDEF,<:DiabaticFrictionCalculator}, t)
    Calculators.update_electronics!(sim.calculator, r)
    for I in eachindex(dv)
        dv[I] = -sim.calculator.adiabatic_derivative[I][1,1]
    end
    DynamicsUtils.divide_by_mass!(dv, sim.atoms.masse)
end

"""
    friction!(g, r, sim, t)

Evaluates friction tensor
"""
function friction!(g, r, sim::AbstractSimulation{<:AbstractMDEF}, t)
    Calculators.evaluate_friction!(sim.calculator, r)
    g .= sim.calculator.friction ./ sim.method.mass_scaling
end

function DynamicsMethods.create_problem(u0, tspan::Tuple, sim::AbstractSimulation{<:AbstractMDEF})
    StochasticDiffEq.DynamicalSDEProblem(acceleration!, DynamicsUtils.velocity!, friction!,
        DynamicsUtils.get_velocities(u0), DynamicsUtils.get_positions(u0), tspan, sim)
end
