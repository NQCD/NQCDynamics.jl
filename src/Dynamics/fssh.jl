export SurfaceHoppingPhasespace
export FSSH
export fssh_callback

struct FSSH{T} <: Method
    nonadiabatic_coupling::Matrix{Matrix{T}}
    density_propagator::Matrix{Complex{T}}
    hopping_probability::Vector{T}
    momentum_rescale::Vector{T}
    function FSSH{T}(DoFs::Integer, atoms::Integer, states::Integer) where {T}
        nonadiabatic_coupling = [zeros(states, states) for i=1:DoFs, j=1:atoms]
        density_propagator = zeros(states, states)
        hopping_probability = zeros(states)
        momentum_rescale = zeros(atoms)
        new{T}(nonadiabatic_coupling, density_propagator, hopping_probability, momentum_rescale)
    end
end

mutable struct SurfaceHoppingPhasespace{T} <: DynamicalVariables{T}
    x::ArrayPartition{Complex{T}, Tuple{Matrix{T}, Matrix{T}, Matrix{Complex{T}}}}
    state::UInt
end

function SurfaceHoppingPhasespace(x::ArrayPartition{T}, state::Integer) where {T<:AbstractFloat}
    SurfaceHoppingPhasespace{T}(ArrayPartition(x.x[1], x.x[2], Complex.(x.x[3])), state)
end

function SurfaceHoppingPhasespace(R::Matrix{T}, P::Matrix{T}, σ::Matrix{Complex{T}}, state::Integer) where {T}
    SurfaceHoppingPhasespace{T}(ArrayPartition(R, P, σ), state)
end

function SurfaceHoppingPhasespace(R::Matrix{T}, P::Matrix{T}, n_states::Integer, state::Integer) where {T}
    σ = zeros(Complex{T}, n_states, n_states)
    σ[state, state] = 1
    SurfaceHoppingPhasespace(R, P, σ, state)
end

get_density_matrix(z::SurfaceHoppingPhasespace) = z.x.x[3]

function create_problem(u0::SurfaceHoppingPhasespace, tspan::Tuple, sim::AbstractSimulation{<:FSSH})
    ODEProblem(motion!, u0, tspan, sim; callback=fssh_callback)
end

function motion!(du::SurfaceHoppingPhasespace, u::SurfaceHoppingPhasespace, sim::Simulation{<:FSSH}, t)
    set_velocity!(du, u, sim)
    update_electronics!(sim, u)
    set_force!(du, u, sim)
    set_density_matrix_derivative!(du, u, sim)
end

function update_electronics!(sim::Simulation{<:FSSH}, u::SurfaceHoppingPhasespace)
    Calculators.evaluate_potential!(sim.calculator, get_positions(u))
    Calculators.evaluate_derivative!(sim.calculator, get_positions(u))
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

function set_force!(du::SurfaceHoppingPhasespace, u::SurfaceHoppingPhasespace, sim::Simulation{<:FSSH})
    for i=1:length(sim.atoms)
        for j=1:sim.DoFs
            get_momenta(du)[j,i] = -sim.calculator.adiabatic_derivative[j,i][u.state, u.state]
        end
    end
#    println(get_momenta(du))
end

function set_density_matrix_derivative!(du::SurfaceHoppingPhasespace, u::SurfaceHoppingPhasespace, sim::Simulation{<:FSSH})
    σ = get_density_matrix(u)
    velocity = get_positions(du)
    V = sim.method.density_propagator

    V .= diagm(sim.calculator.eigenvalues)
    for i=1:length(sim.atoms)
        for j=1:sim.DoFs
            V .-= im*velocity[j,i].*sim.method.nonadiabatic_coupling[j,i]
        end
    end
    get_density_matrix(du) .= -im*(V*σ - σ*V)
end

condition(u, t, integrator::DiffEqBase.DEIntegrator) = true

function affect!(integrator::DiffEqBase.DEIntegrator)
    update_hopping_probability!(integrator)
    
    new_state = select_new_state(integrator.p.method.hopping_probability, integrator.u.state)
    
    if new_state != 0
        if calculate_rescaling_constant!(integrator, new_state)
            execute_hop!(integrator, new_state)
        end
    end
end

function update_hopping_probability!(integrator::DiffEqBase.DEIntegrator)
    sim = integrator.p
    coupling = sim.method.nonadiabatic_coupling
    velocity = get_positions(get_du(integrator))
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

function calculate_rescaling_constant!(integrator::DiffEqBase.DEIntegrator, new_state)::Bool
    sim = integrator.p
    old_state = integrator.u.state
    velocity = get_positions(get_du(integrator))
    
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
        integrator.p.method.momentum_rescale .= min.(abs.((b .+ root) ./ a), abs.((b .- root) ./ a))
        return true
    end
end

function calculate_potential_energy_change(eigenvalues::Vector, new_state::Integer, current_state::Integer)
    eigenvalues[new_state] - eigenvalues[current_state]
end

function execute_hop!(integrator::DiffEqBase.DEIntegrator, new_state::Integer)
    for i in range(integrator.p.atoms)
        coupling = [integrator.p.method.nonadiabatic_coupling[j,i][new_state, integrator.u.state] for j=1:integrator.p.DoFs]
        get_momenta(integrator.u)[:,i] .-= integrator.p.method.momentum_rescale[i] .* coupling
    end
    integrator.u.state = new_state
end

fssh_callback = DiscreteCallback(condition, affect!; save_positions=(false, false))
