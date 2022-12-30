using NQCDistributions: PureState, Diabatic
using NQCDynamics: get_temperature

mutable struct CME{T} <: SurfaceHopping
    hopping_probability::T
    state::Int
    new_state::Int
end

function CME{T}(nstates::Integer) where {T}
    nstates == 2 || error("CME only works with two states at the moment.")
    hopping_probability = zero(T)
    state = 0
    new_state = 0
    return CME(hopping_probability, state, new_state)
end

function Simulation{CME}(atoms::Atoms{T}, model; kwargs...) where {T}
    Simulation(atoms, model, CME{T}(NQCModels.nstates(model)); kwargs...)
end

function DynamicsMethods.motion!(du, u, sim::Simulation{<:CME}, t)
    dr = DynamicsUtils.get_positions(du)
    dv = DynamicsUtils.get_velocities(du)

    r = DynamicsUtils.get_positions(u)
    v = DynamicsUtils.get_velocities(u)

    set_state!(u, sim.method.state) # Make sure the state variables match, 

    DynamicsUtils.velocity!(dr, v, r, sim, t) # Set the velocity
    DynamicsUtils.acceleration!(dv, v, r, sim, t, sim.method.state) # Set the acceleration
end

function DynamicsMethods.DynamicsVariables(::Simulation{<:CME}, v, r, electronic::PureState{Diabatic})
    return SurfaceHoppingVariables(ComponentVector(v=v, r=r), electronic.state)
end

function DynamicsUtils.acceleration!(dv, v, r, sim::Simulation{<:CME}, t, state)
    derivative = Calculators.get_derivative(sim.calculator, r)
    for I in eachindex(dv)
        dv[I] = -derivative[I][state, state]
    end
    DynamicsUtils.divide_by_mass!(dv, masses(sim))
    return nothing
end

function evaluate_hopping_probability!(sim::Simulation{<:CME}, u, dt)
    r = DynamicsUtils.get_positions(u)
    V = Calculators.get_potential(sim.calculator, r)
    ΔV = V[2,2] - V[1,1]
    Γ = 2π * V[2,1]^2
    f = DynamicsUtils.fermi(ΔV, NQCModels.fermilevel(sim), 1/get_temperature(sim))

    if sim.method.state == 1
        sim.method.hopping_probability = Γ * f * dt
    elseif sim.method.state == 2
        sim.method.hopping_probability = Γ * (1 - f) * dt
    end
end
    
function select_new_state(sim::Simulation{<:CME}, u)::Int
    random_number = rand()
    if random_number < sim.method.hopping_probability
        if sim.method.state == 1
            return 2
        elseif sim.method.state == 2
            return 1
        end
    end
    return sim.method.state
end

function rescale_velocity!(::Simulation{<:CME}, u)::Bool
    return true
end

function DynamicsUtils.classical_potential_energy(sim::Simulation{<:CME}, u)
    V = Calculators.get_potential(sim.calculator, DynamicsUtils.get_positions(u))
    return V[u.state, u.state]
end

function DynamicsUtils.classical_hamiltonian(sim::Simulation{<:CME}, u)
    kinetic = DynamicsUtils.classical_kinetic_energy(sim, DynamicsUtils.get_velocities(u))
    potential = DynamicsUtils.classical_potential_energy(sim, u)
    return kinetic + potential
end
