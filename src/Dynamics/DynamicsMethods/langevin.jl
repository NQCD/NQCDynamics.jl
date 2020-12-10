export Langevin
export random_force!

struct Langevin{T<:AbstractFloat} <: Systems.DynamicsParameters
    η::T
end
Langevin(η::Integer) = Langevin{Float64}(η)

function System{Langevin}(atomic_parameters::AtomicParameters{T}, model::Models.Model, 
                temperature::Unitful.Temperature{<:Real}, η::Real, n_DoF::Integer=3) where {T}
    System{Langevin, T}(n_DoF, austrip(temperature), atomic_parameters, model, Langevin(η))
end

function set_force!(du::Phasespace, u::Phasespace, p::System{Langevin})
    Electronics.calculate_derivative!(p.model, p.electronics, get_positions(u))
    get_momenta(du) .= -p.electronics.D0 .- p.dynamics.η .* get_momenta(u)
end

function set_force!(du::RingPolymerPhasespace, u::RingPolymerPhasespace, p::RingPolymerSystem{Langevin})
    for i=1:n_beads(p)
        Electronics.calculate_derivative!(p.model, p.electronics[i], get_positions(u, i))
        get_momenta(du, i) .= -p.electronics[i].D0 .- p.dynamics.η .* get_momenta(u, i)
    end
    apply_interbead_coupling!(du, u, p)
end

function random_force!(du::Phasespace, u::Phasespace, p::System{Langevin}, t)
    get_momenta(du) .= sqrt.(2p.dynamics.η * p.temperature  .* masses(p)')
end
