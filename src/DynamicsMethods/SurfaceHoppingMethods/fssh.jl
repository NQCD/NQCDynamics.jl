using StatsBase: mean
using .NonadiabaticDistributions: NonadiabaticDistribution
using NonadiabaticMolecularDynamics: masses

export FSSH

"""
    FSSH{T} <: SurfaceHopping

Type for fewest-switches surface hopping

```jldoctest
Simulation{FSSH}(Atoms(:H), Free())

# output

Simulation{FSSH{Float64}}:
  Atoms{1, Float64}([:H], UInt8[0x01], [1837.4715941070515])
  Free(1)
```
"""
mutable struct FSSH{T} <: SurfaceHopping
    density_propagator::Matrix{Complex{T}}
    hopping_probability::Vector{T}
    state::Int
    new_state::Int
    rescaling::Symbol
    function FSSH{T}(states::Integer, rescaling::Symbol) where {T}
        density_propagator = zeros(states, states)
        hopping_probability = zeros(states)
        new{T}(density_propagator, hopping_probability, 0, 0, rescaling)
    end
end

function Simulation{FSSH}(atoms::Atoms{S,T}, model::Model; rescaling=:standard, kwargs...) where {S,T}
    Simulation(atoms, model, FSSH{T}(NonadiabaticModels.nstates(model), rescaling); kwargs...)
end

function DynamicsMethods.DynamicsVariables(
    sim::AbstractSimulation{<:SurfaceHopping}, v, r, electronic::NonadiabaticDistribution
)
    σ = NonadiabaticDistributions.initialise_adiabatic_density_matrix(electronic, sim.calculator, r)
    state = sample(Weights(diag(real.(σ))))
    return SurfaceHoppingVariables(ComponentVector(v=v, r=r, σreal=σ, σimag=zero(σ)), state)
end

function acceleration!(dv, v, r, sim::AbstractSimulation{<:FSSH}, t, state)
    for I in eachindex(dv)
        dv[I] = -sim.calculator.adiabatic_derivative[I][state, state]
    end
    DynamicsUtils.divide_by_mass!(dv, masses(sim))
    return nothing
end

"""
    evaluate_hopping_probability!(sim::Simulation{<:FSSH}, u, dt)

Evaluates the probability of hopping from the current state to all other states

# Implementation

- `σ` is Hermitan so the choice `σ[m,s]` or `σ[s,m]` is irrelevant; we take the real part.
- 'd' is skew-symmetric so here the indices are important.
"""
function evaluate_hopping_probability!(sim::AbstractSimulation{<:FSSH}, u, dt)
    v = get_hopping_velocity(sim, u)
    σ = DynamicsUtils.get_quantum_subsystem(u)
    s = sim.method.state
    d = get_hopping_nonadiabatic_coupling(sim)

    fewest_switches_probability!(sim.method.hopping_probability, v, σ, s, d, dt)
end

function get_hopping_nonadiabatic_coupling(sim::Simulation{<:SurfaceHopping})
    sim.calculator.nonadiabatic_coupling
end

function get_hopping_velocity(::Simulation{<:SurfaceHopping}, u)
    DynamicsUtils.get_velocities(u)
end

function get_hopping_eigenvalues(sim::Simulation{<:SurfaceHopping})
    sim.calculator.eigenvalues
end

function fewest_switches_probability!(probability, v, σ, s, d, dt)
    probability .= 0 # Set all entries to 0
    for m in axes(σ, 1)
        if m != s
            for I in eachindex(v)
                probability[m] += 2v[I]*real(σ[m,s]/σ[s,s])*d[I][s,m] * dt
            end
        end
    end

    clamp!(probability, 0, 1) # Restrict probabilities between 0 and 1
    cumsum!(probability, probability)
end

function select_new_state(sim::AbstractSimulation{<:FSSH}, u)

    random_number = rand()
    for (i, prob) in enumerate(sim.method.hopping_probability)
        if i != sim.method.state # Avoid self-hops
            if prob > random_number
                return i
            end
        end
    end
    return sim.method.state
end

"""
    rescale_velocity!(sim::AbstractSimulation{<:FSSH}, u)::Bool

Rescale the velocity in the direction of the nonadiabatic coupling.

# References

[HammesSchiffer1994](@cite)
"""
function rescale_velocity!(sim::AbstractSimulation{<:FSSH}, u)::Bool
    sim.method.rescaling === :off && return true

    new_state, old_state = unpack_states(sim)
    velocity = get_hopping_velocity(sim, u)
    eigs = get_hopping_eigenvalues(sim)
    
    d = extract_nonadiabatic_coupling(get_hopping_nonadiabatic_coupling(sim), new_state, old_state)
    a = calculate_a(sim, d)
    b = calculate_b(d, velocity)
    c = calculate_potential_energy_change(eigs, new_state, old_state)

    discriminant = b^2 - 4a * c
    discriminant < 0 && return false # Frustrated hop with insufficient kinetic energy

    root = sqrt(discriminant)
    if b < 0
        γ = (b + root) / 2a
    else
        γ = (b - root) / 2a
    end
    perform_rescaling!(sim, DynamicsUtils.get_velocities(u), γ, d)

    return true
end

"""
    unpack_states(sim)

Get the two states that we are hopping between.
"""
function unpack_states(sim::AbstractSimulation{<:FSSH})
    return (sim.method.new_state, sim.method.state)
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
        a += coupling[I]^2 * masses(sim, I)
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

function calculate_potential_energy_change(eigs::AbstractVector, new_state::Integer, current_state::Integer)
    return eigs[new_state] - eigs[current_state]
end

function Estimators.diabatic_population(sim::AbstractSimulation{<:FSSH}, u)
    r = DynamicsUtils.get_positions(u)
    U = NonadiabaticDistributions.evaluate_transformation(sim.calculator, r)

    σ = copy(DynamicsUtils.get_quantum_subsystem(u).re)
    σ[diagind(σ)] .= 0
    σ[u.state, u.state] = 1

    return diag(U * σ * U')
end

function Estimators.adiabatic_population(sim::AbstractSimulation{<:FSSH}, u)
    population = zeros(NonadiabaticModels.nstates(sim.calculator.model))
    population[u.state] = 1
    return population
end

function DynamicsUtils.classical_hamiltonian(sim::Simulation{<:FSSH}, u)
    kinetic = DynamicsUtils.classical_kinetic_energy(sim, DynamicsUtils.get_velocities(u))
    Calculators.evaluate_potential!(sim.calculator, DynamicsUtils.get_positions(u))
    Calculators.eigen!(sim.calculator)
    potential = sim.calculator.eigenvalues[sim.method.state]
    return kinetic + potential
end
