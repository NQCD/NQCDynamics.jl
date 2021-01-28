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

function set_force!(du::Phasespace, u::Phasespace, sim::Simulation{<:AbstractMDEF})
    Calculators.evaluate_derivative!(sim.calculator, get_positions(u))
    get_momenta(du) .= -sim.calculator.derivative
    
    Calculators.evaluate_friction!(sim.calculator, get_positions(u))
    mul!(sim.method.drag, sim.calculator.friction, get_flat_momenta(u))

    drag = reshape(sim.method.drag, sim.DoFs, length(sim.atoms))
    
    get_momenta(du) .-= drag ./ sim.atoms.masses'
end

function random_force!(du, u::Phasespace, sim::Simulation{<:AbstractMDEF}, t)
    slice = sim.DoFs*length(sim.atoms)+1
    du[slice:end, slice:end] .= sqrt.(sim.calculator.friction .* 2 .* get_temperature(sim, t))
end

get_temperature(sim::Simulation{<:MDEF}, ::AbstractFloat) = sim.temperature
get_temperature(sim::Simulation{<:TwoTemperatureMDEF}, t::AbstractFloat) = sim.method.temperature(t)

function create_problem(u0::Phasespace, tspan::Tuple, sim::Simulation{<:AbstractMDEF})
    n = sim.DoFs * length(sim.atoms) * 2
    SDEProblem(motion!, random_force!, u0, tspan, sim; noise_rate_prototype=zeros(n,n))
end
select_algorithm(::AbstractSimulation{<:MDEF}) = SRA()
