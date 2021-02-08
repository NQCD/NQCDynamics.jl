using Random
using LinearAlgebra

export MDEF

abstract type AbstractMDEF <: Method end

struct MDEF{T<:AbstractFloat} <: AbstractMDEF
    drag::Vector{T}
    function MDEF{T}(atoms::Integer; DoF::Integer=3) where {T}
        new(zeros(DoF*atoms))
    end
end

struct TwoTemperatureMDEF{T<:AbstractFloat} <: AbstractMDEF
    drag::Vector{T}
    temperature::Function
    function TwoTemperatureMDEF{T}(atoms::Integer, temperature::Function; DoF::Integer=3) where {T}
        new(zeros(DoF*atoms), temperature)
    end
end

get_temperature(sim::Simulation{<:MDEF}, ::AbstractFloat) = sim.temperature
get_temperature(sim::Simulation{<:TwoTemperatureMDEF}, t::AbstractFloat) = sim.method.temperature(t)

function friction!(du, r, sim, t)
    Calculators.evaluate_friction!(sim.calculator, r)

    du.x[1] .= sim.calculator.friction
    du.x[2] .= sqrt.(get_temperature(sim, t) ./ repeat(sim.atoms.masses; inner=sim.DoFs))
end

function create_problem(u0::ClassicalDynamicals, tspan::Tuple, sim::AbstractSimulation{<:AbstractMDEF})
    DynamicalSDEProblem(acceleration!, velocity!, friction!, get_velocities(u0), get_positions(u0), tspan, sim)
end

# function create_problem(u0::ArrayPartition, tspan::Tuple, sim::AbstractSimulation{<:Classical})
#     DynamicalODEProblem(acceleration!, velocity!, u0.x[1], u0.x[2], tspan, sim)
# end

select_algorithm(::AbstractSimulation{<:AbstractMDEF}) = MDEF_BAOAB()
