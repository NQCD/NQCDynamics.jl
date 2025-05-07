using NQCDistributions: PureState, Diabatic
using NQCDynamics: get_temperature
using QuadGK: QuadGK

abstract type ClassicalMasterEquation <: SurfaceHopping end

function DynamicsMethods.motion!(du, u, sim::Simulation{<:ClassicalMasterEquation}, t)
    dr = DynamicsUtils.get_positions(du)
    dv = DynamicsUtils.get_velocities(du)

    r = DynamicsUtils.get_positions(u)
    v = DynamicsUtils.get_velocities(u)

    set_state!(u, sim.method.state) # Make sure the state variables match, 

    DynamicsUtils.velocity!(dr, v, r, sim, t) # Set the velocity
    DynamicsUtils.acceleration!(dv, v, r, sim, t) # Set the acceleration
end

function DynamicsMethods.DynamicsVariables(::AbstractSimulation{<:ClassicalMasterEquation}, v, r, electronic::PureState{Diabatic})
    nt = (x = ComponentVector(v=v, r=r), state = electronic.state)
    return NamedArrayPartition(nt)
end

function evaluate_hopping_probability!(sim::Simulation{<:ClassicalMasterEquation}, u, dt)
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
    
function select_new_state(sim::AbstractSimulation{<:ClassicalMasterEquation}, u)::Int
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

function rescale_velocity!(::AbstractSimulation{<:ClassicalMasterEquation}, u)::Bool
    return true
end

"""
    CME{T} <: ClassicalMasterEquation

Simple surface hopping method for Newns-Anderson (Anderson-Holstein) models.

- Dou, Nitzan, Subotnik, J. Chem. Phys. 142, 084110 (2015)
- Dou, Subotnik, J. Phys. Chem. A, 24, 757-771 (2020)
"""
mutable struct CME{T} <: ClassicalMasterEquation
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

function DynamicsUtils.acceleration!(dv, v, r, sim::AbstractSimulation{<:CME}, t)
    derivative = Calculators.get_derivative(sim.calculator, r)
    state = sim.method.state
    for I in eachindex(dv, derivative)
        dv[I] = -derivative[I][state, state]
    end
    DynamicsUtils.divide_by_mass!(dv, masses(sim))
    return nothing
end

function DynamicsUtils.classical_potential_energy(sim::Simulation{<:CME}, u)
    V = Calculators.get_potential(sim.calculator, DynamicsUtils.get_positions(u))
    int_state = convert(Int, u.state)
    return V[int_state, int_state]
end

"""
    BCME{T} <: ClassicalMasterEquation

Extension to CME that incorporates broadening in the potential energy surfaces.

Note that we do not rescale the velocity as this is not mentioned only in the 2016 paper,
not any of the later ones, so I presume they later decided not to do it and instead keep it the same as the original CME.

- Dou, Subotnik, J. Chem. Phys. 144, 024116 (2016)
- Dou, Subotnik, J. Phys. Chem. A, 24, 757-771 (2020)
"""
mutable struct BCME{T} <: ClassicalMasterEquation
    hopping_probability::T
    state::Int
    new_state::Int
    bandwidth::T
end

function BCME{T}(nstates::Integer, bandwidth) where {T}
    nstates == 2 || error("BCME only works with two states at the moment.")
    hopping_probability = zero(T)
    state = 0
    new_state = 0
    return BCME(hopping_probability, state, new_state, bandwidth)
end

function Simulation{BCME}(atoms::Atoms{T}, model; bandwidth, kwargs...) where {T}
    Simulation(atoms, model, BCME{T}(NQCModels.nstates(model), bandwidth); kwargs...)
end

function DynamicsUtils.acceleration!(dv, v, r, sim::Simulation{<:BCME}, t)
    state = sim.method.state
    ∂V = Calculators.get_derivative(sim.calculator, r)
    V = Calculators.get_potential(sim.calculator, r)
    Γ = 2π * V[2,1]^2
    β = 1 / get_temperature(sim, t)
    μ = NQCModels.fermilevel(sim)
    h = V[2,2] - V[1,1]
    n = evaluate_broadening(h, μ, β, Γ, sim.method.bandwidth)
    n2 = evaluate_force_broadening(h, μ, β, Γ, sim.method.bandwidth)
    f = DynamicsUtils.fermi(h, μ, β)

    for I in eachindex(dv, ∂V)
        ∂U0 = ∂V[I][1,1]
        ∂Γᵢ = 4π * V[2,1] * ∂V[I][2,1]
        ∂U1 = ∂V[I][2,2]
        ∂h = ∂U1 - ∂U0
        F₁ = -∂h*n
        F₂ = -∂Γᵢ/Γ*n2
        dv[I] = - ∂V[I][state,state] + f*∂h + F₁ + F₂
    end
    DynamicsUtils.divide_by_mass!(dv, masses(sim))
    return nothing
end

"""
    evaluate_broadening(h, μ, β, Γ)

Evaluate the convolution of the Fermi function with a Lorentzian.
"""
function evaluate_broadening(h, μ, β, Γ, bandwidth)
    A(ϵ) = 1/π * Γ/2 / ((ϵ-h)^2 + (Γ/2)^2)
    kernel(ϵ) = A(ϵ) * DynamicsUtils.fermi(ϵ, μ, β)

    integral, _ = QuadGK.quadgk(kernel, μ-bandwidth/2, μ+bandwidth/2)

    return integral
end

function evaluate_force_broadening(h, μ, β, Γ, bandwidth)
    A(ϵ) = 1/π * Γ/2 / ((ϵ-h)^2 + (Γ/2)^2)
    kernel(ϵ) = (ϵ - h) * A(ϵ) * DynamicsUtils.fermi(ϵ, μ, β)

    integral, _ = QuadGK.quadgk(kernel, μ-bandwidth/2, μ+bandwidth/2)

    return integral
end
