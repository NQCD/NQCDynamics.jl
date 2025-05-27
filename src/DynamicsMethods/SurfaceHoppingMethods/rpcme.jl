
function DynamicsMethods.motion!(du, u, sim::RingPolymerSimulation{<:ClassicalMasterEquation}, t)
    dr = DynamicsUtils.get_positions(du)
    dv = DynamicsUtils.get_velocities(du)

    r = DynamicsUtils.get_positions(u)
    v = DynamicsUtils.get_velocities(u)

    set_state!(u, sim.method.state) # Make sure the state variables match, 

    DynamicsUtils.velocity!(dr, v, r, sim, t) # Set the velocity
    DynamicsUtils.acceleration!(dv, v, r, sim, t) # Set the acceleration
    DynamicsUtils.apply_interbead_coupling!(dv, r, sim)
end

function evaluate_hopping_probability!(sim::RingPolymerSimulation{<:ClassicalMasterEquation}, u, dt)
    r = DynamicsUtils.get_positions(u)
    V = Calculators.get_centroid_potential(sim.calculator, r)
    ΔV = V[2,2] - V[1,1]
    Γ = 2π * V[2,1]^2
    f = DynamicsUtils.fermi(ΔV, NQCModels.fermilevel(sim), 1/get_temperature(sim))

    if sim.method.state == 1
        sim.method.hopping_probability = Γ * f * dt
    elseif sim.method.state == 2
        sim.method.hopping_probability = Γ * (1 - f) * dt
    end
end

function RingPolymerSimulation{CME}(atoms::Atoms{T}, model, n_beads; kwargs...) where {T}
    RingPolymerSimulation(atoms, model, CME{T}(NQCModels.nstates(model)), n_beads; kwargs...)
end

function DynamicsUtils.classical_potential_energy(sim::RingPolymerSimulation{<:CME}, u)
    r = DynamicsUtils.get_positions(u)
    V = Calculators.get_potential(sim.calculator, r)

    # pre-accessing state from u - seems to fix a Int rounding error with rpcme phonon relaxation test 
    state = u.state

    potential = zero(eltype(r))
    @debug state_rounding_error = first(state) - round(Int, first(state))
    int_state = round(Int, first(state))
    for b in axes(r,3) # eachbead
        potential += V[b][int_state, int_state]
    end
    return potential
end

function RingPolymerSimulation{BCME}(atoms::Atoms{T}, model, n_beads; bandwidth, kwargs...) where {T}
    RingPolymerSimulation(atoms, model, BCME{T}(NQCModels.nstates(model), bandwidth), n_beads; kwargs...)
end

function DynamicsUtils.acceleration!(dv, v, r, sim::RingPolymerSimulation{<:BCME}, t)
    ∂V = Calculators.get_derivative(sim.calculator, r)
    V = Calculators.get_potential(sim.calculator, r)
    state = sim.method.state

    β = 1 / get_temperature(sim, t)
    μ = NQCModels.fermilevel(sim)

    for b in axes(r,3) # eachbead
        Γ = 2π * V[b][2,1]^2
        h = V[b][2,2] - V[b][1,1]
        n = evaluate_broadening(h, μ, β, Γ, sim.method.bandwidth)
        n2 = evaluate_force_broadening(h, μ, β, Γ, sim.method.bandwidth)
        f = DynamicsUtils.fermi(h, μ, β)

        for i in NQCModels.mobileatoms(sim)
            for j in dofs(sim)
                ∂U0 = ∂V[j,i,b][1,1]
                ∂Γᵢ = 4π * V[b][2,1] * ∂V[j,i,b][2,1]
                ∂U1 = ∂V[j,i,b][2,2]
                ∂h = ∂U1 - ∂U0
                F₁ = -∂h*n
                F₂ = -∂Γᵢ/Γ*n2
                dv[j,i,b] = - ∂V[j,i,b][state,state] + f*∂h + F₁ + F₂
            end
        end
    end

    DynamicsUtils.divide_by_mass!(dv, masses(sim))
    return nothing
end
