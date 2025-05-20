
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

function NQCDynamics.Simulation{MDEF}(atoms::Atoms, model::Model; kwargs...)
    NQCDynamics.Simulation(atoms, model, MDEF(atoms.masses, ndofs(model));  kwargs...)
end
function NQCDynamics.RingPolymerSimulation{MDEF}(atoms::Atoms, model::Model, n_beads::Integer; kwargs...)
    NQCDynamics.RingPolymerSimulation(atoms, model, MDEF(atoms.masses, ndofs(model)), n_beads; kwargs...)
end

"""
    friction!(g, r, sim, t)

Evaluates friction tensor
"""
function friction!(g, r, sim::AbstractSimulation{<:AbstractMDEF}, t)
    friction = NQCCalculators.get_friction(sim.cache, r)
    g .= friction ./ sim.method.mass_scaling
end

function DynamicsMethods.create_problem(u0, tspan::Tuple, sim::AbstractSimulation{<:AbstractMDEF})
    StochasticDiffEq.DynamicalSDEProblem(acceleration!, DynamicsUtils.velocity!, friction!,
        DynamicsUtils.get_velocities(u0), DynamicsUtils.get_positions(u0), tspan, sim)
end
