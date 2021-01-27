export SurfaceHoppingPhasespace
export IESH
export IESH_callback
using Revise

struct IESH{T} <: Method
    nonadiabatic_coupling::Matrix{Matrix{T}}
    density_propagator::Matrix{Complex{T}}
    hopping_probability::Vector{T}
    momentum_rescale::Vector{T}
    function IESH{T}(DoFs::Integer, atoms::Integer, states::Integer) where {T}
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

# construct the density matrix
function SurfaceHoppingPhasespace(R::Matrix{T}, P::Matrix{T}, n_states::Integer, state::Integer) where {T}
    σ = zeros(Complex{T}, n_states, n_states)
    σ[state, state] = 1
#    println(σ)
    SurfaceHoppingPhasespace(R, P, σ, state)
end

get_density_matrix(z::SurfaceHoppingPhasespace) = z.x.x[3]

# This belongs to run_trajectory.
# It first follows motion! and then iesh_callback
function create_problem(u0::SurfaceHoppingPhasespace, tspan::Tuple, sim::AbstractSimulation{<:IESH})
    ODEProblem(motion!, u0, tspan, sim; callback=IESH_callback)
end

function motion!(du::SurfaceHoppingPhasespace, u::SurfaceHoppingPhasespace, sim::Simulation{<:IESH}, t)
    set_velocity!(du, u, sim)  # classical.jl momenta/(atom_mass) v = p/m 
    update_electronics!(sim, u) # The next three routines 
    set_force!(du, u, sim) # 4th routine
    set_density_matrix_derivative!(du, u, sim)# 5th routine
end

# Get the current eneries, the current derivate and the current eigenvectors
# see ../Calculators/Calculators.jl for the functions used here.
# Then get the eigenvectors of the energy
# then use them to transform the derivatives into the adiabatic form.
function update_electronics!(sim::Simulation{<:IESH}, u::SurfaceHoppingPhasespace)
    Calculators.evaluate_potential!(sim.calculator, get_positions(u))
    Calculators.evaluate_derivative!(sim.calculator, get_positions(u))
    Calculators.eigen!(sim.calculator)
    Calculators.transform_derivative!(sim.calculator)
    evaluate_nonadiabatic_coupling!(sim)
end

#These routines point to the previous routine and hands down an empty array,
# the adiabatic_derivatives calculated in transform_drivatives and the 
# energy eigenvectors
function evaluate_nonadiabatic_coupling!(sim::Simulation{<:IESH})
    evaluate_nonadiabatic_coupling!.(
        sim.method.nonadiabatic_coupling,
        sim.calculator.adiabatic_derivative,
        Ref(sim.calculator.eigenvalues))
end

# Calculates the nonadiabatic coupling. This is probably also true for IESH
function evaluate_nonadiabatic_coupling!(coupling::Matrix, adiabatic_derivative::Matrix, eigenvalues::Vector)
    for i=1:length(eigenvalues)
        for j=i+1:length(eigenvalues)
            coupling[j,i] = -adiabatic_derivative[j,i] / (eigenvalues[j]-eigenvalues[i])
            coupling[i,j] = -coupling[j,i]
        end
    end
end

# gets the momenta corresponding to the current state from the adiabatic derivates 
# (transformation of the diabatic ones, see 3 routines above)
function set_force!(du::SurfaceHoppingPhasespace, u::SurfaceHoppingPhasespace, sim::Simulation{<:IESH})
    for i=1:length(sim.atoms)
        for j=1:sim.DoFs
            get_momenta(du)[j,i] = -sim.calculator.adiabatic_derivative[j,i][u.state, u.state]
        end
    end
end

# Calculate the _time_-derivate  of the  density matrix
# Goes back to motion!
function set_density_matrix_derivative!(du::SurfaceHoppingPhasespace, u::SurfaceHoppingPhasespace, sim::Simulation{<:IESH})
    σ = get_density_matrix(u)
    velocity = get_positions(du)
    V = sim.method.density_propagator

    V .= diagm(sim.calculator.eigenvalues)# Creates a diagonal matrix from eigenvalues
    # Calculation going on here from: Martens_JPhysChemA_123_1110_2019, eq. 6
    # d is the nonadiabatic coupling matrix
    #i ħ dσ/dt = iħ sum_l [(V_{m,l} - i ħ velocity d_{m,l})*σ_{l,n} - &
    #                      σ_{m,l}*(V_{l,n} - iħ velocity d_{l,n})]
    # vgl. also SubotnikBellonzi_AnnuRevPhyschem_67_387_2016, eq. 5 and 6
    for i=1:length(sim.atoms)
        for j=1:sim.DoFs
            V .-= im*velocity[j,i].*sim.method.nonadiabatic_coupling[j,i]
        end
    end
    get_density_matrix(du) .= -im*(V*σ - σ*V)
end

# This sets when the condition will be true. In this case, always.
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
                    println(velocity[j,i]," ", real(σ[m,s]/σ[s,s])," ",coupling[j,i][s,m]," ", dt)
                end
            end
        end
    end
    println(clamp!(sim.method.hopping_probability, 0, 1))
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

IESH_callback = DiscreteCallback(condition, affect!; save_positions=(false, false))
