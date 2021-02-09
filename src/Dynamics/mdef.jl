using Random
using LinearAlgebra

export MDEF
export TwoTemperatureMDEF

abstract type AbstractMDEF <: Method end

struct MDEF <: AbstractMDEF end

struct TwoTemperatureMDEF <: AbstractMDEF
    temperature::Function
end

get_temperature(sim::Simulation{MDEF}, ::AbstractFloat) = sim.temperature
get_temperature(sim::Simulation{TwoTemperatureMDEF}, t::AbstractFloat) = sim.method.temperature(t)

function friction!(du, r, sim, t)
    Calculators.evaluate_friction!(sim.calculator, r)

    du.x[1] .= sim.calculator.friction
    du.x[2] .= sqrt.(get_temperature(sim, t) ./ repeat(sim.atoms.masses; inner=sim.DoFs))
end

function create_problem(u0::ClassicalDynamicals, tspan::Tuple, sim::AbstractSimulation{<:AbstractMDEF})
    DynamicalSDEProblem(acceleration!, velocity!, friction!, get_velocities(u0), get_positions(u0), tspan, sim)
end

select_algorithm(::AbstractSimulation{<:AbstractMDEF}) = MDEF_BAOAB()
