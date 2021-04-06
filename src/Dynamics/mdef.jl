using Random
using LinearAlgebra
using ..Calculators: DiabaticFrictionCalculator

export MDEF
export TwoTemperatureMDEF

abstract type AbstractMDEF <: Method end

"""
$(TYPEDEF)

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

"""
$(TYPEDEF)

Same as standard MDEF but uses a function to determine the time-dependent temperature.
"""
struct TwoTemperatureMDEF{M} <: AbstractMDEF
    mass_scaling::M
    temperature::Function
end

TwoTemperatureMDEF(masses::AbstractVector, DoFs::Integer, temperature::Function) =
    TwoTemperatureMDEF(get_mass_scale_matrix(masses, DoFs), temperature)

function get_mass_scale_matrix(masses::AbstractVector, DoFs::Integer)
    vect = repeat(sqrt.(masses), inner=DoFs)
    vect * vect'
end

"""Gets the temperature as a function of time during MDEF."""
get_temperature(sim::Simulation{<:MDEF}, ::AbstractFloat) = sim.temperature
get_temperature(sim::Simulation{<:TwoTemperatureMDEF}, t::AbstractFloat) = sim.method.temperature(t)

"""
    acceleration!(dv, v, r, sim::Simulation{MDEF,<:DiabaticFrictionCalculator}, t)

Sets acceleration due to ground state force when using a `DiabaticFrictionModel`.
"""
function acceleration!(dv, v, r, sim::Simulation{MDEF,<:DiabaticFrictionCalculator}, t)
    Calculators.update_electronics!(sim.calculator, r)
    for i in axes(r, 2)
        for j in axes(r, 1)
            dv[j,i] = -sim.calculator.adiabatic_derivative[j,i][1,1] / sim.atoms.masses[i]
        end
    end
end

"""
    friction!(du, r, sim, t)

Evaluates friction tensor and provides variance of random force.
"""
function friction!(du, r, sim::AbstractSimulation{<:AbstractMDEF}, t)
    Calculators.evaluate_friction!(sim.calculator, r)

    du.x[1] .= sim.calculator.friction ./ sim.method.mass_scaling
    du.x[2] .= sqrt.(get_temperature(sim, t) ./ repeat(sim.atoms.masses; inner=sim.DoFs))
end

function create_problem(u0::ClassicalDynamicals, tspan::Tuple, sim::AbstractSimulation{<:AbstractMDEF})
    create_problem(u0.x, tspan, sim)
end

function create_problem(u0::ArrayPartition, tspan::Tuple, sim::AbstractSimulation{<:AbstractMDEF})
    DynamicalSDEProblem(acceleration!, velocity!, friction!, get_velocities(u0), get_positions(u0), tspan, sim)
end

select_algorithm(::AbstractSimulation{<:AbstractMDEF}) = MDEF_BAOAB()
