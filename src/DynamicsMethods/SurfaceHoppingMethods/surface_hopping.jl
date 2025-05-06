
function check_hop!(u, t, integrator)::Bool
    sim = integrator.p
    evaluate_hopping_probability!(sim, u, OrdinaryDiffEq.get_proposed_dt(integrator))
    set_new_state!(sim.method, select_new_state(sim, u))
    return sim.method.new_state != sim.method.state
end

function execute_hop!(integrator)
    sim = integrator.p
    if rescale_velocity!(sim, integrator.u)
        set_state!(integrator.u, sim.method.new_state)
        set_state!(sim.method, sim.method.new_state)
    end
    return nothing
end

set_state!(container, new_state::Integer) = container.state = new_state
set_state!(container, new_state::AbstractVector) = copyto!(container.state, new_state)
set_new_state!(container, new_state::Integer) = container.new_state = new_state
set_new_state!(container, new_state::AbstractVector) = copyto!(container.new_state, new_state)

const HoppingCallback = DiffEqBase.DiscreteCallback(check_hop!, execute_hop!;
                                                    save_positions=(false, false))

DynamicsMethods.get_callbacks(::AbstractSimulation{<:SurfaceHopping}) = HoppingCallback

"""
This function should set the field `sim.method.hopping_probability`.
"""
function evaluate_hopping_probability!(::AbstractSimulation{<:SurfaceHopping}, u, dt)
    throw(error("Implement this for your method."))
end

"""
This function should return the desired state determined by the probability.
Should return the original state if no hop is desired.
"""
function select_new_state(::AbstractSimulation{<:SurfaceHopping}, u)
    throw(error("Implement this for your method."))
end

"""
    rescale_velocity!(sim::AbstractSimulation{<:SurfaceHopping}, u)::Bool

Rescale the velocity in the direction of the nonadiabatic coupling.

# References

[HammesSchiffer1994](@cite)
"""
function rescale_velocity!(sim::AbstractSimulation{<:SurfaceHopping}, u)::Bool
    sim.method.rescaling === :off && return true #no rescaling so always accept hop

    new_state, old_state = unpack_states(sim)
    velocity = DynamicsUtils.get_hopping_velocity(sim, u.v)
    r = Du.r
    eigs = DynamicsUtils.get_hopping_eigenvalues(sim, r)
    
    d = extract_nonadiabatic_coupling(DynamicsUtils.get_hopping_nonadiabatic_coupling(sim, r), new_state, old_state)
    a = calculate_a(sim, d)
    b = calculate_b(d, velocity)
    c = calculate_potential_energy_change(eigs, new_state, old_state)

    discriminant = b^2 - 4a * c
    if discriminant < 0 # frustrated hop
        if sim.method.rescaling == :standard
            # Frustrated hop with insufficient kinetic energy, no velocity inversion
            # Follow 1990 Tully recipe and do nothing
            return false 
        elseif sim.method.rescaling == :vinversion
            # Frustrated hop with insufficient kinetic energy
            # perform inversion of the velocity component along the nonadiabatic coupling vector
            frustrated_hop_invert_velocity!(sim, u.v, d)
            return false
        else   
            throw(error("This mode of rescaling is not implemented"))
        end
    else     #sufficient energy for hopping
        root = sqrt(discriminant)
        if b < 0
            γ = (b + root) / 2a
        else
            γ = (b - root) / 2a
        end
        perform_rescaling!(sim, u.v, γ, d)

        return true
    end
end

"""
    extract_nonadiabatic_coupling(coupling, new_state, old_state)

Extract the nonadiabatic coupling vector between states `new_state` and `old_state`
"""
function extract_nonadiabatic_coupling(coupling, new_state, old_state)
    [coupling[I][new_state, old_state] for I ∈ CartesianIndices(coupling)]
end

"""
    calculate_a(sim::AbstractSimulation{<:SurfaceHopping}, coupling::AbstractMatrix)

Equation 40 from [HammesSchiffer1994](@cite).
"""
function calculate_a(sim::AbstractSimulation{<:SurfaceHopping}, coupling::AbstractMatrix)
    a = zero(eltype(coupling))
    for I in CartesianIndices(coupling)
        a += coupling[I]^2 / masses(sim, I)
    end
    return a / 2
end

"""
    calculate_b(coupling::AbstractMatrix, velocity::AbstractMatrix)

Equation 41 from [HammesSchiffer1994](@cite).
"""
function calculate_b(coupling::AbstractMatrix, velocity::AbstractMatrix)
    return dot(coupling, velocity)
end

"""
    perform_rescaling!(
        sim::AbstractSimulation{<:SurfaceHopping}, velocity, velocity_rescale, d
    )

Equation 33 from [HammesSchiffer1994](@cite).
"""
function perform_rescaling!(
    sim::AbstractSimulation{<:SurfaceHopping}, velocity, γ, d
)
    for I in CartesianIndices(d)
        velocity[I] -= γ * d[I] / masses(sim, I)
    end
    return nothing
end

"""
    frustrated_hop_invert_velocity!(
        sim::AbstractSimulation{<:SurfaceHopping}, velocity, d
    )

Measures the component of velocity along the nonadiabatic coupling vector and inverts that component.
"""
function frustrated_hop_invert_velocity!(
    sim::AbstractSimulation{<:SurfaceHopping}, velocity, d
)
    dn = LinearAlgebra.normalize(d)
    γ = dot(velocity,dn)
    for I in CartesianIndices(dn)
        velocity[I] -= 2γ * dn[I]
    end
    return nothing
end

function calculate_potential_energy_change(eigs::AbstractVector, new_state::Integer, current_state::Integer)
    return eigs[new_state] - eigs[current_state]
end
