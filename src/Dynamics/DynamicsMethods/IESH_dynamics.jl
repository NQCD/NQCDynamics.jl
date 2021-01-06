export IESH_dynamics

# This file lists the un-changable parts for the IESH-dynamics and 
# if it makes sense???
# some functions for calling specific information like forces or moment. 

struct IESH_dynamics{T<:AbstractFloat} <: Systems.DynamicsParameters
    #η::T
end
IESH_dynamics() = IESH_dynamics{Float64}()

function System{IESH_dynamics}(atomic_parameters::AtomicParameters{T}, model::Models.Model, 
    temperature::Unitful.Temperature{<:Real}, n_DoF::Integer=3) where {T}
System{IESH_dynamics, T}(n_DoF, austrip(temperature), atomic_parameters, model, IESH_dynamics())
end

# Not sure if this is necessary or if it can be called directly form the model
#function set_force!(du::Phasespace, u::Phasespace, p::System{Langevin})
#    Electronics.calculate_derivative!(p.model, p.electronics, get_positions(u))
#    get_momenta(du) .= -p.electronics.D0 .- p.dynamics.η .* get_momenta(u)
#end





