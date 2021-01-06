using Random
using LinearAlgebra

export MDEF

struct MDEF{T<:AbstractFloat} <: Method
    drag::Vector{T}
    function MDEF{T}(atoms::Integer, DoF::Integer) where {T}
        new(zeros(DoF*atoms))
    end
end

function set_force!(du::Phasespace, u::Phasespace, sim::Simulation{<:MDEF})
    Calculators.evaluate_derivative!(sim.calculator, get_positions(u))
    get_momenta(du) .= -sim.calculator.derivative

    Calculators.evaluate_friction!(sim.calculator, get_positions(u))
    mul!(sim.method.drag, sim.calculator.friction, get_flat_momenta(u))
    
    get_flat_momenta(du) .-= sim.method.drag
end

function random_force!(du, u::Phasespace, sim::Simulation{<:MDEF}, t)
    Calculators.evaluate_friction!(sim.calculator, get_positions(u))
    du[sim.DoFs*length(sim.atoms)+1:end, sim.DoFs*length(sim.atoms)+1:end] .= sim.calculator.friction * sqrt(2sim.temperature)
end

function create_problem(u0::Phasespace, tspan::Tuple, sim::Simulation{<:MDEF})
    n = sim.DoFs * length(sim.atoms) * 2
    SDEProblem(motion!, random_force!, u0, tspan, sim; noise_rate_prototype=zeros(n,n))
end
select_algorithm(::AbstractSimulation{<:MDEF}) = LambaEM()
