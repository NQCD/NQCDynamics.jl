export Langevin
export random_force!

struct Langevin{T<:AbstractFloat} <: Systems.DynamicsParameters
    temperature::T
    η::T
end

Langevin(temperature::Real, η::Real) = Langevin{Float64}(temperature, η)

function Langevin(temperature::Unitful.Temperature{<:Real}, η::Real)
    Langevin(austrip(temperature), η)
end

function System{Langevin}(atomic_parameters::AtomicParameters, model::Models.Model, 
                temperature::Unitful.Temperature{<:Real}, η::Real, n_DoF::Integer=3)
    System{Langevin}(n_DoF, atomic_parameters, model, Langevin(temperature, η))
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
    get_momenta(du) .= sqrt.(2p.dynamics.η * p.dynamics.temperature  .* masses(p)')
end
