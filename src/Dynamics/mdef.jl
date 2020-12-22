using Random
using LinearAlgebra

export MDEF

struct MDEF{T<:AbstractFloat} <: Method
    friction::Matrix{T}
    drag::Vector{T}
    function MDEF{T}(atoms::Integer, DoF::Integer) where {T}
        new(zeros(T, atoms*DoF, atoms*DoF), zeros(DoF*atoms))
    end
end

function set_force!(du::Phasespace, u::Phasespace, sim::Simulation{<:MDEF})
    get_momenta(du) .= -Calculators.evaluate_derivative(sim.calculator, get_positions(u))

    update_friction!(sim.method.friction)
    mul!(sim.method.drag, sim.method.friction, get_flat_momenta(u))
    
    get_flat_momenta(du) .-= sim.method.drag
end

function random_force!(du, u::Phasespace, sim::Simulation{<:MDEF}, t)
    update_friction!(sim.method.friction)
    du[sim.DoFs*length(sim.atoms)+1:end, sim.DoFs*length(sim.atoms)+1:end] .= sim.method.friction * sqrt(2sim.temperature)
end

update_friction!(friction) = randn!(friction)
