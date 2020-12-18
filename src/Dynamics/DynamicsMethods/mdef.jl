using Random
using LinearAlgebra

export MDEF

struct MDEF{T<:AbstractFloat} <: Systems.DynamicsParameters
    friction::Matrix{T}
    drag::Vector{T}
    function MDEF{T}(atoms::Integer, DoF::Integer) where {T}
        new(zeros(T, atoms*DoF, atoms*DoF), zeros(DoF*atoms))
    end
end

function System{MDEF}(atomic_parameters::AtomicParameters{T}, model::Models.Model, 
                temperature::Unitful.Temperature{<:Real}, n_DoF::Integer=3) where {T}
    System{MDEF,T}(n_DoF, austrip(temperature), atomic_parameters, model,
        MDEF{T}(atomic_parameters.n_atoms, n_DoF))
end

function set_force!(du::Phasespace, u::Phasespace, p::System{MDEF})
    Electronics.calculate_derivative!(p.model, p.electronics, get_positions(u))
    get_momenta(du) .= -p.electronics.D0

    update_friction!(p.dynamics.friction)
    mul!(p.dynamics.drag, p.dynamics.friction, get_flat_momenta(u))
    
    get_flat_momenta(du) .-= p.dynamics.drag
end

function random_force!(du, u::Phasespace, p::System{MDEF}, t)
    update_friction!(p.dynamics.friction)
    du[n_DoF(p)*n_atoms(p)+1:end, n_DoF(p)*n_atoms(p)+1:end] .= p.dynamics.friction
end

update_friction!(friction) = randn!(friction)
