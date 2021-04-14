using StatsBase: mean
using .Calculators: DiabaticCalculator, RingPolymerDiabaticCalculator

export SurfaceHoppingDynamicals
export FSSH
export fssh_callback
export get_density_matrix

abstract type SurfaceHopping <: Method end

struct FSSH{T} <: SurfaceHopping
    density_propagator::Matrix{Complex{T}}
    hopping_probability::Vector{T}
    new_state::Vector{Int}
    function FSSH{T}(states::Integer) where {T}
        density_propagator = zeros(states, states)
        hopping_probability = zeros(states)
        new_state = [0]
        new{T}(density_propagator, hopping_probability, new_state)
    end
end

mutable struct SurfaceHoppingDynamicals{T,D}  <: DynamicalVariables{T}
    x::ArrayPartition{Complex{T}, Tuple{D,D,Matrix{Complex{T}}}}
    state::Int
end

function SurfaceHoppingDynamicals(x::ArrayPartition{T}, state::Integer) where {T<:AbstractFloat}
    SurfaceHoppingDynamicals(ArrayPartition(x.x[1], x.x[2], Complex.(x.x[3])), state)
end

function SurfaceHoppingDynamicals(v::AbstractArray, r::AbstractArray, σ::Matrix{Complex{T}}, state::Integer) where {T}
    SurfaceHoppingDynamicals(ArrayPartition(v, r, σ), state)
end

function SurfaceHoppingDynamicals(v::AbstractArray, r::AbstractArray, n_states::Integer, state::Integer)
    σ = zeros(Complex{eltype(r)}, n_states, n_states)
    σ[state, state] = 1
    SurfaceHoppingDynamicals(v, r, σ, state)
end

get_density_matrix(u::SurfaceHoppingDynamicals) = u.x.x[3]

function create_problem(u0::SurfaceHoppingDynamicals, tspan::Tuple, sim::AbstractSimulation{<:FSSH})
    ODEProblem(motion!, u0, tspan, sim)
end

get_callbacks(::AbstractSimulation{<:FSSH}) = fssh_callback

function motion!(du::DynamicalVariables, u::DynamicalVariables, sim::AbstractSimulation{<:SurfaceHopping}, t)
    dr = get_positions(du)
    dv = get_velocities(du)
    dσ = get_density_matrix(du)

    r = get_positions(u)
    v = get_velocities(u)
    σ = get_density_matrix(u)

    velocity!(dr, v, r, sim, t)
    update_electronics!(sim, r)
    acceleration!(dv, v, r, sim, t, u.state)
    set_density_matrix_derivative!(dσ, v, σ, sim)
end

function update_electronics!(sim::AbstractSimulation{<:SurfaceHopping}, r::AbstractArray)
    Calculators.evaluate_potential!(sim.calculator, r)
    Calculators.evaluate_derivative!(sim.calculator, r)
    Calculators.eigen!(sim.calculator)
    Calculators.transform_derivative!(sim.calculator)
    Calculators.evaluate_nonadiabatic_coupling!(sim.calculator)
end

function acceleration!(dv, v, r, sim::Simulation{<:FSSH}, t, state=1)
    for i=1:length(sim.atoms)
        for j=1:sim.DoFs
            dv[j,i] = -sim.calculator.adiabatic_derivative[j,i][state, state] / sim.atoms.masses[i]
        end
    end
end

function acceleration!(dv, v, r, sim::RingPolymerSimulation{<:FSSH}, t, state=1)
    for i in axes(v, 3)
        for j=1:length(sim.atoms)
            for k=1:sim.DoFs
                dv[k,j,i] = -sim.calculator.adiabatic_derivative[k,j,i][state, state] / sim.atoms.masses[j]
            end
        end
    end
    apply_interbead_coupling!(dv, r, sim)
end

function set_density_matrix_derivative!(dσ, v, σ, sim::Simulation{<:SurfaceHopping})
    V = sim.method.density_propagator

    V .= diagm(sim.calculator.eigenvalues)
    for i=1:length(sim.atoms)
        for j=1:sim.DoFs
            V .-= im*v[j,i].*sim.calculator.nonadiabatic_coupling[j,i]
        end
    end
    dσ .= -im*(V*σ - σ*V)
end

function set_density_matrix_derivative!(dσ, v, σ, sim::RingPolymerSimulation{<:SurfaceHopping})
    V = sim.method.density_propagator

    V .= diagm(mean(sim.calculator.eigenvalues))
    ivd = mean(im .* v .* sim.calculator.nonadiabatic_coupling; dims=3)
    for I in eachindex(ivd)
        V .-= ivd[I]
    end
    dσ .= -im*(V*σ - σ*V)
end

function condition!(u, t, integrator)::Bool
    update_hopping_probability!(integrator.p, integrator.u, get_proposed_dt(integrator))
    select_new_state!(integrator.p, integrator.u)
    integrator.p.method.new_state[1] != 0
end

function execute_hop!(integrator)
    rescale_velocity!(integrator) && (integrator.u.state = integrator.p.method.new_state[1])
end

function update_hopping_probability!(sim::AbstractSimulation{<:FSSH}, u, dt)
    coupling = sim.calculator.nonadiabatic_coupling
    velocity = get_velocities(u)
    s = u.state
    σ = get_density_matrix(u)
    
    update_hopping_probability!(sim, s, σ, velocity, coupling, dt)

    clamp!(sim.method.hopping_probability, 0, 1) # Restrict probabilities between 0 and 1
    cumsum!(sim.method.hopping_probability, sim.method.hopping_probability)
end

function update_hopping_probability!(sim::Simulation{<:FSSH}, s, σ, velocity, coupling, dt)
    sim.method.hopping_probability .= 0 # Set all entries to 0
    for m=1:sim.calculator.model.n_states
        if m != s
            for I in eachindex(velocity)
                sim.method.hopping_probability[m] += 2velocity[I]*real(σ[m,s]/σ[s,s])*coupling[I][s,m] * dt
            end
        end
    end
end

function update_hopping_probability!(sim::RingPolymerSimulation{<:FSSH}, s, σ, velocity, coupling, dt)
    sim.method.hopping_probability .= 0 # Set all entries to 0
    for m=1:sim.calculator.model.n_states
        if m != s
            for I in eachindex(velocity)
                sim.method.hopping_probability[m] += 2velocity[I]*real(σ[m,s]/σ[s,s])*coupling[I][s,m] * dt
            end
        end
    end
    sim.method.hopping_probability ./= length(sim.beads)
end

function select_new_state!(sim::AbstractSimulation{<:FSSH}, u)
    sim.method.new_state[1] = select_new_state(sim.method.hopping_probability, u.state)
end

function select_new_state(probability::Vector{T}, current_state::Integer)::UInt where {T<:AbstractFloat}
    random_number = rand()
    for (i, prob) in enumerate(probability)
        if i != current_state # Avoid self-hops
            if prob > random_number
                return i # Return index of selected state
            end
        end
    end
    0 # Return 0 if no hop is desired
end

function rescale_velocity!(integrator::DiffEqBase.DEIntegrator)::Bool
    sim = integrator.p
    old_state = integrator.u.state
    new_state = integrator.p.method.new_state[1]
    velocity = get_velocities(integrator.u)
    
    c = calculate_potential_energy_change(sim.calculator, new_state, old_state)
    a, b = evaluate_a_and_b(sim, velocity, new_state, old_state)
    discriminant = b.^2 .- 2a.*c

    any(discriminant .< 0) && return false

    root = sqrt.(discriminant)
    velocity_rescale = min.(abs.((b .+ root) ./ a), abs.((b .- root) ./ a))
    perform_rescaling!(sim, velocity, velocity_rescale, new_state, old_state)

    return true
end

function evaluate_a_and_b(sim::AbstractSimulation{<:FSSH}, velocity::AbstractArray, new_state, old_state)
    a = zeros(length(sim.atoms))
    b = zero(a)
    @views for i in range(sim.atoms)
        coupling = [sim.calculator.nonadiabatic_coupling[j,i][new_state, old_state] for j=1:sim.DoFs]
        a[i] = coupling'coupling / sim.atoms.masses[i]
        b[i] = velocity[:,i]'coupling
    end
    (a, b)
end

function evaluate_a_and_b(sim::RingPolymerSimulation{<:FSSH}, velocity::RingPolymerArray, new_state, old_state)
    a = zeros(length(sim.atoms), length(sim.beads))
    b = zero(a)
    @views for i in range(sim.beads)
        for j in range(sim.atoms)
            coupling = [sim.calculator.nonadiabatic_coupling[k,j,i][new_state, old_state] for k=1:sim.DoFs]
            a[j,i] = coupling'coupling / sim.atoms.masses[j]
            b[j,i] = velocity[:,j,i]'coupling
        end
    end
    (a, b)
end

function perform_rescaling!(sim::Simulation{<:FSSH}, velocity, velocity_rescale, new_state, old_state)
    for i in range(sim.atoms)
        coupling = [sim.calculator.nonadiabatic_coupling[j,i][new_state, old_state] for j=1:sim.DoFs]
        velocity[:,i] .-= velocity_rescale[i] .* coupling ./ sim.atoms.masses[i]
    end
end

function perform_rescaling!(sim::RingPolymerSimulation{<:FSSH}, velocity, velocity_rescale, new_state, old_state)
    for i=1:length(sim.beads)
        for j=1:length(sim.atoms)
            coupling = [sim.calculator.nonadiabatic_coupling[k,j,i][new_state, old_state] for k=1:sim.DoFs]
            velocity[:,j,i] .-= velocity_rescale[j,i] .* coupling ./ sim.atoms.masses[j]
        end
    end
end

function calculate_potential_energy_change(calc::DiabaticCalculator, new_state::Integer, current_state::Integer)
    calc.eigenvalues[new_state] - calc.eigenvalues[current_state]
end

function calculate_potential_energy_change(calc::RingPolymerDiabaticCalculator, new_state::Integer, current_state::Integer)
    mean([eigs[new_state] - eigs[current_state] for eigs in calc.eigenvalues])
end

const fssh_callback = DiscreteCallback(condition!, execute_hop!; save_positions=(false, false))
