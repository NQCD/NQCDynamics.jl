using StatsBase: mean
using .Calculators: AbstractDiabaticCalculator

export FSSH

mutable struct FSSH{T} <: SurfaceHopping
    density_propagator::Matrix{Complex{T}}
    hopping_probability::Vector{T}
    state::Int
    new_state::Int
    function FSSH{T}(states::Integer) where {T}
        density_propagator = zeros(states, states)
        hopping_probability = zeros(states)
        new{T}(density_propagator, hopping_probability, 0, 0)
    end
end

function Dynamics.DynamicsVariables(sim::Simulation{<:SurfaceHopping}, v, r, state::Integer; type=:diabatic)
    n_states = sim.calculator.model.n_states
    if type == :diabatic
        Calculators.evaluate_potential!(sim.calculator, r)
        Calculators.eigen!(sim.calculator)
        U = sim.calculator.eigenvectors

        diabatic_density = zeros(n_states, n_states)
        diabatic_density[state, state] = 1
        σ = U' * diabatic_density * U
        state = sample(Weights(diag(real.(σ))))
    else
        σ = zeros(n_states, n_states)
        σ[state, state] = 1
    end
    return SurfaceHoppingVariables(ComponentVector(v=v, r=r, σreal=σ, σimag=zero(σ)), state)
end

function acceleration!(dv, v, r, sim::Simulation{<:FSSH}, t, state)
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
    v = Dynamics.get_velocities(u)
    σ = Dynamics.get_quantum_subsystem(u)
    s = sim.method.state
    d = sim.calculator.nonadiabatic_coupling

    fewest_switches_probability!(sim.method.hopping_probability, v, σ, s, d, dt)
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
    old_state = sim.method.state
    new_state = sim.method.new_state
    velocity = Dynamics.get_velocities(u)
    
    c = calculate_potential_energy_change(sim.calculator, new_state, old_state)
    a, b = evaluate_a_and_b(sim, velocity, new_state, old_state)
    discriminant = b.^2 .- 2a.*c

    any(discriminant .< 0) && return false

    root = sqrt.(discriminant)
    plus = (b .+ root) ./ a
    minus = (b .- root) ./ a 
    velocity_rescale = sum(abs.(plus)) < sum(abs.(minus)) ? plus : minus
    perform_rescaling!(sim, velocity, velocity_rescale, new_state, old_state)

    return true
end

function evaluate_a_and_b(sim::Simulation{<:SurfaceHopping}, velocity, new_state, old_state)
    a = zeros(length(sim.atoms))
    b = zero(a)
    @views for i in range(sim.atoms)
        coupling = [sim.calculator.nonadiabatic_coupling[j,i][new_state, old_state] for j=1:sim.DoFs]
        a[i] = dot(coupling, coupling) / sim.atoms.masses[i]
        b[i] = dot(velocity[:,i], coupling)
    end
    return (a, b)
end

function perform_rescaling!(sim::Simulation{<:SurfaceHopping}, velocity, velocity_rescale, new_state, old_state)
    for i in range(sim.atoms)
        coupling = [sim.calculator.nonadiabatic_coupling[j,i][new_state, old_state] for j=1:sim.DoFs]
        velocity[:,i] .-= velocity_rescale[i] .* coupling ./ sim.atoms.masses[i]
    end
    return nothing
end

function calculate_potential_energy_change(calc::AbstractDiabaticCalculator, new_state::Integer, current_state::Integer)
    return calc.eigenvalues[new_state] - calc.eigenvalues[current_state]
end

"""
    get_diabatic_population(sim::Simulation{<:FSSH}, u)

Diabatic populations recommended from J. Chem. Phys. 139, 211101 (2013).
"""
function Estimators.diabatic_population(sim::Simulation{<:FSSH}, u)
    Calculators.evaluate_potential!(sim.calculator, get_positions(u))
    U = eigvecs(sim.calculator.potential)

    σ = copy(Dynamics.get_quantum_subsystem(u).re)
    σ[diagind(σ)] .= 0
    σ[u.state, u.state] = 1

    return diag(U * σ * U')
end

"""
    get_adiabatic_population(sim::AbstractSimulation{<:FSSH}, u)

Adiabatic population directly from discrete state variable.
"""
function Estimators.adiabatic_population(sim::AbstractSimulation{<:FSSH}, u)
    population = zeros(sim.calculator.model.n_states)
    population[u.state] = 1
    return population
end

function NonadiabaticMolecularDynamics.evaluate_hamiltonian(sim::Simulation{<:FSSH}, u)
    k = NonadiabaticMolecularDynamics.evaluate_kinetic_energy(sim.atoms.masses, get_velocities(u))
    Calculators.evaluate_potential!(sim.calculator, get_positions(u))
    Calculators.eigen!(sim.calculator)
    p = sim.calculator.eigenvalues[sim.method.state]
    return k + p
end