using Random
using LinearAlgebra
using ..Calculators: DiabaticFrictionCalculator

export MDEF
export TwoTemperatureMDEF

abstract type AbstractMDEF <: Method end

struct MDEF <: AbstractMDEF end

struct TwoTemperatureMDEF <: AbstractMDEF
    temperature::Function
end

"""Gets the temperature as a function of time during MDEF."""
get_temperature(sim::Simulation{MDEF}, ::AbstractFloat) = sim.temperature
get_temperature(sim::Simulation{TwoTemperatureMDEF}, t::AbstractFloat) = sim.method.temperature(t)

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
function friction!(du, r, sim, t)
    Calculators.evaluate_friction!(sim.calculator, r)

    du.x[1] .= sim.calculator.friction
    if sim isa Simulation{MDEF,<:DiabaticFrictionCalculator}
        du.x[1] ./= sim.atoms.masses[1]
    end
    du.x[2] .= sqrt.(get_temperature(sim, t) ./ repeat(sim.atoms.masses; inner=sim.DoFs))
end

function create_problem(u0::ClassicalDynamicals, tspan::Tuple, sim::AbstractSimulation{<:AbstractMDEF})
    create_problem(u0.x, tspan, sim)
end

function create_problem(u0::ArrayPartition, tspan::Tuple, sim::AbstractSimulation{<:AbstractMDEF})
    DynamicalSDEProblem(acceleration!, velocity!, friction!, get_velocities(u0), get_positions(u0), tspan, sim)
end

select_algorithm(::AbstractSimulation{<:AbstractMDEF}) = MDEF_BAOAB()
