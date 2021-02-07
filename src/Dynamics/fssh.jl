export SurfaceHoppingDynamicals
export FSSH
export fssh_callback
export get_density_matrix

struct FSSH{T} <: Method
    nonadiabatic_coupling::Matrix{Matrix{T}}
    density_propagator::Matrix{Complex{T}}
    hopping_probability::Vector{T}
    velocity_rescale::Vector{T}
    new_state::Vector{Int}
    function FSSH{T}(DoFs::Integer, atoms::Integer, states::Integer) where {T}
        nonadiabatic_coupling = [zeros(states, states) for i=1:DoFs, j=1:atoms]
        density_propagator = zeros(states, states)
        hopping_probability = zeros(states)
        velocity_rescale = zeros(atoms)
        new_state = [0]
        new{T}(nonadiabatic_coupling, density_propagator, hopping_probability, velocity_rescale, new_state)
    end
end

mutable struct SurfaceHoppingDynamicals{T}  <: DynamicalVariables{T}
    x::ArrayPartition{Complex{T}, Tuple{Matrix{T}, Matrix{T}, Matrix{Complex{T}}}}
    state::UInt
end

function SurfaceHoppingDynamicals(x::ArrayPartition{T}, state::Integer) where {T<:AbstractFloat}
    SurfaceHoppingDynamicals{T}(ArrayPartition(x.x[1], x.x[2], Complex.(x.x[3])), state)
end

function SurfaceHoppingDynamicals(v::Matrix{T}, r::Matrix{T}, σ::Matrix{Complex{T}}, state::Integer) where {T}
    SurfaceHoppingDynamicals{T}(ArrayPartition(v, r, σ), state)
end

function SurfaceHoppingDynamicals(v::Matrix{T}, r::Matrix{T}, n_states::Integer, state::Integer) where {T}
    σ = zeros(Complex{T}, n_states, n_states)
    σ[state, state] = 1
    SurfaceHoppingDynamicals(v, r, σ, state)
end

get_density_matrix(u::SurfaceHoppingDynamicals) = u.x.x[3]

function create_problem(u0::SurfaceHoppingDynamicals, tspan::Tuple, sim::AbstractSimulation{<:FSSH})
    ODEProblem(motion!, u0, tspan, sim; callback=fssh_callback)
end

function motion!(du::SurfaceHoppingDynamicals, u::SurfaceHoppingDynamicals, sim::Simulation{<:FSSH}, t)
    dr = get_positions(du)
    dv = get_velocities(du)
    dσ = get_density_matrix(du)

    r = get_positions(u)
    v = get_velocities(u)
    σ = get_density_matrix(u)

    velocity!(dr, v, r, sim, t)
    update_electronics!(sim, r)
    acceleration!(dv, v, r, sim, t; state=u.state)
    set_density_matrix_derivative!(dσ, v, σ, sim)
end

function update_electronics!(sim::Simulation{<:FSSH}, r::Matrix)
    Calculators.evaluate_potential!(sim.calculator, r)
    Calculators.evaluate_derivative!(sim.calculator, r)
    Calculators.eigen!(sim.calculator)
    Calculators.transform_derivative!(sim.calculator)
    evaluate_nonadiabatic_coupling!(sim)
end

function evaluate_nonadiabatic_coupling!(sim::Simulation{<:FSSH})
    evaluate_nonadiabatic_coupling!.(
        sim.method.nonadiabatic_coupling,
        sim.calculator.adiabatic_derivative,
        Ref(sim.calculator.eigenvalues))
end

function evaluate_nonadiabatic_coupling!(coupling::Matrix, adiabatic_derivative::Matrix, eigenvalues::Vector)
    for i=1:length(eigenvalues)
        for j=i+1:length(eigenvalues)
            coupling[j,i] = -adiabatic_derivative[j,i] / (eigenvalues[j]-eigenvalues[i])
            coupling[i,j] = -coupling[j,i]
        end
    end
end

function acceleration!(dv, v, r, sim::Simulation{<:FSSH}, t; state=1)
    for i=1:length(sim.atoms)
        for j=1:sim.DoFs
            dv[j,i] = -sim.calculator.adiabatic_derivative[j,i][state, state] / sim.atoms.masses[i]
        end
    end
end

function set_density_matrix_derivative!(dσ, v, σ, sim::Simulation{<:FSSH})
    V = sim.method.density_propagator

    V .= diagm(sim.calculator.eigenvalues)
    for i=1:length(sim.atoms)
        for j=1:sim.DoFs
            V .-= im*v[j,i].*sim.method.nonadiabatic_coupling[j,i]
        end
    end
    dσ .= -im*(V*σ - σ*V)
end

function condition!(u, t, integrator)
    update_hopping_probability!(integrator)
    select_new_state!(integrator)
    integrator.p.method.new_state[1] != 0
end

function affect!(integrator::DiffEqBase.DEIntegrator)
    calculate_rescaling_constant!(integrator) && execute_hop!(integrator)
end

function update_hopping_probability!(integrator::DiffEqBase.DEIntegrator)
    sim = integrator.p
    coupling = sim.method.nonadiabatic_coupling
    velocity = get_velocities(integrator.u)
    s = integrator.u.state
    σ = get_density_matrix(integrator.u)
    dt = get_proposed_dt(integrator)
    
    sim.method.hopping_probability .= 0 # Set all entries to 0
    for m=1:integrator.p.calculator.model.n_states
        if m != s
            for i=1:length(sim.atoms)
                for j=1:sim.DoFs
                    sim.method.hopping_probability[m] += 2velocity[j,i]*real(σ[m,s]/σ[s,s])*coupling[j,i][s,m] * dt
                end
            end
        end
    end
    clamp!(sim.method.hopping_probability, 0, 1) # Restrict probabilities between 0 and 1
    cumsum!(sim.method.hopping_probability, sim.method.hopping_probability)
end

function select_new_state!(integrator)
    integrator.p.method.new_state[1] = select_new_state(integrator.p.method.hopping_probability, integrator.u.state)
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

function calculate_rescaling_constant!(integrator::DiffEqBase.DEIntegrator)::Bool
    sim = integrator.p
    old_state = integrator.u.state
    new_state = integrator.p.method.new_state[1]
    velocity = get_velocities(integrator.u)
    
    c = calculate_potential_energy_change(sim.calculator.eigenvalues, new_state, old_state)
    
    a = zeros(length(sim.atoms))
    b = zero(a)
    @views for i in range(sim.atoms)
        coupling = [sim.method.nonadiabatic_coupling[j,i][new_state, old_state] for j=1:sim.DoFs]
        a[i] = coupling'coupling / sim.atoms.masses[i]
        b[i] = velocity[:,i]'coupling
    end
    
    discriminant = b.^2 .- 2a.*c
    if any(discriminant .< 0)
        return false
    else
        root = sqrt.(discriminant)
        integrator.p.method.velocity_rescale .= min.(abs.((b .+ root) ./ a), abs.((b .- root) ./ a))
        return true
    end
end

function calculate_potential_energy_change(eigenvalues::Vector, new_state::Integer, current_state::Integer)
    eigenvalues[new_state] - eigenvalues[current_state]
end

function execute_hop!(integrator::DiffEqBase.DEIntegrator)
    new_state = integrator.p.method.new_state[1]
    for i in range(integrator.p.atoms)
        coupling = [integrator.p.method.nonadiabatic_coupling[j,i][new_state, integrator.u.state] for j=1:integrator.p.DoFs]
        get_velocities(integrator.u)[:,i] .-= integrator.p.method.velocity_rescale[i] .* coupling ./ integrator.p.atoms.masses[i]
    end
    integrator.u.state = new_state
end

const fssh_callback = DiscreteCallback(condition!, affect!; save_positions=(false, false))
