using Random
using LinearAlgebra
using ..Calculators: DiabaticFrictionCalculator

export MDEF
export TwoTemperatureMDEF

abstract type AbstractMDEF <: Method end

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

"""
    acceleration!(dv, v, r, sim::Simulation{MDEF,<:DiabaticFrictionCalculator}, t)

Sets acceleration due to ground state force when using a `DiabaticFrictionModel`.
"""
function acceleration!(dv, v, r, sim::Simulation{MDEF,<:DiabaticFrictionCalculator}, t)
    Calculators.update_electronics!(sim.calculator, r)
    for I in eachindex(dv)
        dv[I] = -sim.calculator.adiabatic_derivative[I][1,1]
    end
    divide_by_mass!(dv, sim.atoms.masse)
end

"""
    friction!(g, r, sim, t)

Evaluates friction tensor
"""
function friction!(g, r, sim::AbstractSimulation{<:AbstractMDEF}, t)
    Calculators.evaluate_friction!(sim.calculator, r)
    g .= sim.calculator.friction ./ sim.method.mass_scaling
end

function create_problem(u0, tspan::Tuple, sim::AbstractSimulation{<:AbstractMDEF})
    DynamicalSDEProblem(acceleration!, velocity!, friction!, get_velocities(u0), get_positions(u0), tspan, sim)
end

select_algorithm(::Simulation{<:AbstractMDEF}) = MDEF_BAOAB()
select_algorithm(::RingPolymerSimulation{<:AbstractMDEF}) = BCOCB()
