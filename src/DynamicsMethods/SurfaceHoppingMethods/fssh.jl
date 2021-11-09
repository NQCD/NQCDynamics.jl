using StatsBase: mean
using .NonadiabaticDistributions: NonadiabaticDistribution

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
    DynamicsUtils.divide_by_mass!(dv, sim.atoms.masses)
    return nothing
end

"""
    evaluate_hopping_probability!(sim::Simulation{<:FSSH}, u, dt)

Evaluates the probability of hopping from the current state to all other states

# Implementation

- `σ` is Hermitan so the choice `σ[m,s]` or `σ[s,m]` is irrelevant; we take the real part.
- 'd' is skew-symmetric so here the indices are important.
"""
function evaluate_hopping_probability!(sim::Simulation{<:FSSH}, u, dt)
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

function rescale_velocity!(sim::AbstractSimulation{<:FSSH}, u)::Bool
    sim.method.rescaling === :off && return true

    old_state = sim.method.state
    new_state = sim.method.new_state
    velocity = get_hopping_velocity(sim, u)
    
    c = calculate_potential_energy_change(sim, new_state, old_state)
    a, b = evaluate_a_and_b(sim, velocity, new_state, old_state)
    discriminant = b.^2 .- 2a.*c

    any(discriminant .< 0) && return false

    root = sqrt.(discriminant)
    plus = (b .+ root) ./ a
    minus = (b .- root) ./ a 
    velocity_rescale = sum(abs.(plus)) < sum(abs.(minus)) ? plus : minus
    perform_rescaling!(sim, DynamicsUtils.get_velocities(u), velocity_rescale, new_state, old_state)

    return true
end

function evaluate_a_and_b(sim::AbstractSimulation{<:SurfaceHopping}, velocity, new_state, old_state)
    a = zeros(length(sim.atoms))
    b = zero(a)
    d = get_hopping_nonadiabatic_coupling(sim)
    @views for i in range(sim.atoms)
        coupling = [d[j,i][new_state, old_state] for j=1:ndofs(sim)]
        a[i] = dot(coupling, coupling) / sim.atoms.masses[i]
        b[i] = dot(velocity[:,i], coupling)
    end
    return (a, b)
end

function perform_rescaling!(sim::AbstractSimulation{<:SurfaceHopping}, velocity, velocity_rescale, new_state, old_state)
    d = get_hopping_nonadiabatic_coupling(sim)
    for i in range(sim.atoms)
        coupling = [d[j,i][new_state, old_state] for j=1:ndofs(sim)]
        velocity[:,i] .-= velocity_rescale[i] .* coupling ./ sim.atoms.masses[i]
    end
    return nothing
end

function calculate_potential_energy_change(sim::AbstractSimulation{<:FSSH}, new_state::Integer, current_state::Integer)
    eigs = get_hopping_eigenvalues(sim)
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
