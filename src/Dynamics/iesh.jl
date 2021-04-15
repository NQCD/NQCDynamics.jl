export IESHPhasespace
export IESH
export IESH_callback
<<<<<<< HEAD
using Revise
using DiffEqCallbacks
using Combinatorics
using TickTock
=======
using DiffEqCallbacks
>>>>>>> c0f6d97171f30e24cabdbf5157e3fd907556a7a1
using LinearAlgebra

"""This module controles how IESH is executed. For a description of IESH, see e.g.
Roy, Shenvi, Tully, J. Chem. Phys. 130, 174716 (2009) and 
Shenvi, Roy,  Tully,J. Chem. Phys. 130, 174107 (2009).
It first needs to be initialized, e.g.
dynam = Dynamics.IESH{Float64}(DoFs, atoms, n_states)
After setting up the simulation (e.g. sim = Simulation(atoms, model, dynam; DoFs=1)),
for a single trajectory, the Phasespace is set up inside this module:
z = SurfaceHoppingPhasespace(r,p, n_states+2, k)  
and the trajectory can be run with:
solution = Dynamics.run_trajectory(z, (0.0, 15000.0), sim)
The latter calls back to DifferentialEquations, but the callback of DifferentialEquations
with the hopping between surfaces is handled here.
"""

<<<<<<< HEAD
struct IESH{T} <: Method
    nonadiabatic_coupling::Matrix{Matrix{T}}
    density_propagator::Matrix{Complex{T}}
    hopping_probability::Vector{T}
    momentum_rescale::Vector{T}
    tmp_matrix1::Matrix{Complex{T}}
    tmp_matrix2::Matrix{Complex{T}}
    tmp_matrix3::Matrix{Complex{T}}
    # matrix as temporary arrays to allocate only once
    function IESH{T}(DoFs::Integer, atoms::Integer, states::Integer) where {T}
        nonadiabatic_coupling = [zeros(states, states) for i=1:DoFs, j=1:atoms]
        density_propagator = zeros(states, states)
        hopping_probability = zeros(3)
        momentum_rescale = zeros(atoms)
        tmp_matrix1 = zeros(states, states)
        tmp_matrix2 = zeros(states, states)
        tmp_matrix3 = zeros(states, states)
        new{T}(nonadiabatic_coupling, density_propagator, hopping_probability, 
               momentum_rescale, tmp_matrix1, tmp_matrix2, tmp_matrix3)
=======
struct IESH{T} <: SurfaceHopping
    density_propagator::Matrix{Complex{T}}
    hopping_probability::Vector{T}
    momentum_rescale::Vector{T}
    function IESH{T}(DoFs::Integer, atoms::Integer, states::Integer) where {T}
        density_propagator = zeros(states, states)
        #hopping_probability = zeros(states)
        hopping_probability = zeros(3)
        momentum_rescale = zeros(atoms)
        new{T}(density_propagator, hopping_probability, momentum_rescale)
>>>>>>> c0f6d97171f30e24cabdbf5157e3fd907556a7a1
    end
end

mutable struct IESHPhasespace{T} <: DynamicalVariables{T}
    x::ArrayPartition{Complex{T}, Tuple{Matrix{T}, Matrix{T}, Matrix{Complex{T}}}}
    state::Vector{Int}
end

function IESHPhasespace(x::ArrayPartition{T}, state::Vector{Int}) where {T<:AbstractFloat}
    IESHPhasespace{T}(ArrayPartition(x.x[1], x.x[2], Complex.(x.x[3])), state)
end

function IESHPhasespace(R::Matrix{T}, P::Matrix{T}, σ::Matrix{Complex{T}}, state::Vector{Int}) where {T}
    IESHPhasespace{T}(ArrayPartition(R, P, σ), state)
end

# construct the density matrix
<<<<<<< HEAD
# Adapted from James's definition from FSSH. We now need states as the Vector |k>
=======
# Adapted from James's definition. We now need states as the Vector |k>
>>>>>>> c0f6d97171f30e24cabdbf5157e3fd907556a7a1
# See: ShenviRoyTully_JChemPhys_130_174107_2009
function IESHPhasespace(R::Matrix{T}, P::Matrix{T}, n_states::Integer, state::Vector{Int}) where {T}
    σ = zeros(Complex{T}, n_states, n_states)
    # 1 needs to be replaced by the info from state, once I've gotten things to run.
    c = 0
    for i=1:length(state)/2
        #c = c + 1
        σ[Int(i), Int(i)] = 1
    end
    IESHPhasespace(R, P, σ, state)
end

get_density_matrix(z::IESHPhasespace) = z.x.x[3]

# This belongs to run_trajectory.
# It first follows motion! and then iesh_callback
<<<<<<< HEAD
=======
# Doesn't work, yet. HERE!
>>>>>>> c0f6d97171f30e24cabdbf5157e3fd907556a7a1
function create_problem(u0::IESHPhasespace, tspan::Tuple, sim::AbstractSimulation{<:IESH})
    #IESH_callback=create_energy_saving_callback(u0,integrator::DiffEqBase.DEIntegrator)
    # IESH_callback=create_saving_callback()
    #cb2 = DiscreteCallback(condition, affect!; save_positions=(false, false))
    cb2 = DiscreteCallback(condition, affect!; save_positions=(false, false))
    #ODEProblem(motion!, u0, tspan, sim; callback=cb2, dt = 1, adaptive = false)
<<<<<<< HEAD
    ODEProblem(motion!, u0, tspan, sim; callback=cb2, save_everystep=false,dense=false, dt = 1, adaptive = false)
    #ODEProblem(motion!, u0, tspan, sim; callback=cb2, save_everystep=false,dense=false)
 end

function motion!(du::IESHPhasespace, u::IESHPhasespace, sim::Simulation{<:IESH}, t)
    set_velocity!(du, u, sim)  # classical.jl momenta/(atom_mass) v = p/m 
    update_electronics!(sim, u) # The next three routines 
    set_force!(du, u, sim) # 4th routine
    set_density_matrix_derivative!(du, u, sim)# 5th routine
end

# Get the current eneries, the current derivate and the current eigenvectors
# see ../Calculators/Calculators.jl for the functions used here.
# Then get the eigenvectors of the energy
# then use them to transform the derivatives into the adiabatic form.
function update_electronics!(sim::Simulation{<:IESH}, u::IESHPhasespace)
    Calculators.evaluate_potential!(sim.calculator, get_positions(u))
    Calculators.evaluate_derivative!(sim.calculator, get_positions(u))
    Calculators.eigen!(sim.calculator)
    #Calculators.transform_derivative!(sim.calculator)
    # Speeds up by a factor of 10
    Calculators.transform_derivative!(sim.calculator, sim.method.tmp_matrix1, sim.method.tmp_matrix2, sim.method.tmp_matrix3)
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
    #println("ping1")
    #println(eigenvalues[1])
    for i=1:length(eigenvalues)
        for j=i+1:length(eigenvalues)
            coupling[j,i] = -adiabatic_derivative[j,i] / (eigenvalues[j]-eigenvalues[i])
            coupling[i,j] = -coupling[j,i]
        end
    end
    #println(coupling)
end
=======
    #ODEProblem(motion!, u0, tspan, sim; callback=cb2, save_everystep=false,dense=false, dt = 10, adaptive = false)
    ODEProblem(motion!, u0, tspan, sim; callback=cb2, save_everystep=false,dense=false)
 end

# function motion!(du::IESHPhasespace, u::IESHPhasespace, sim::Simulation{<:IESH}, t)
#     #println(get_positions(u))
#     set_velocity!(du, u, sim)  # classical.jl momenta/(atom_mass) v = p/m 
#     update_electronics!(sim, u) # The next three routines 
#     set_force!(du, u, sim) # 4th routine
#     set_density_matrix_derivative!(du, u, sim)# 5th routine
# end

# # Get the current eneries, the current derivate and the current eigenvectors
# # see ../Calculators/Calculators.jl for the functions used here.
# # Then get the eigenvectors of the energy
# # then use them to transform the derivatives into the adiabatic form.
# function update_electronics!(sim::Simulation{<:IESH}, u::IESHPhasespace)
#     Calculators.evaluate_potential!(sim.calculator, get_positions(u))
#     Calculators.evaluate_derivative!(sim.calculator, get_positions(u))
#     Calculators.eigen!(sim.calculator)
#     Calculators.transform_derivative!(sim.calculator)
#     evaluate_nonadiabatic_coupling!(sim)
# end

# #These routines point to the previous routine and hands down an empty array,
# # the adiabatic_derivatives calculated in transform_drivatives and the 
# # energy eigenvectors
# function evaluate_nonadiabatic_coupling!(sim::Simulation{<:IESH})
#     evaluate_nonadiabatic_coupling!.(
#         sim.method.nonadiabatic_coupling,
#         sim.calculator.adiabatic_derivative,
#         Ref(sim.calculator.eigenvalues))
# end

# # Calculates the nonadiabatic coupling. This is probably also true for IESH
# function evaluate_nonadiabatic_coupling!(coupling::Matrix, adiabatic_derivative::Matrix, eigenvalues::Vector)
#     for i=1:length(eigenvalues)
#         for j=i+1:length(eigenvalues)
#             coupling[j,i] = -adiabatic_derivative[j,i] / (eigenvalues[j]-eigenvalues[i])
#             coupling[i,j] = -coupling[j,i]
#         end
#     end
# end
>>>>>>> c0f6d97171f30e24cabdbf5157e3fd907556a7a1

# gets the momenta corresponding to the current state from the adiabatic derivates 
# (transformation of the diabatic ones, see 3 routines above)
# Here, instead of just calculating the moment, I assume I will need to some up (?)
# the different contributions to the force according to Eq.(12) in the Tully paper.
<<<<<<< HEAD
function set_force!(du::IESHPhasespace, u::IESHPhasespace, sim::Simulation{<:IESH})
    for i=1:length(sim.atoms)
        for j=1:sim.DoFs
            get_momenta(du)[j,i] = 0.0
            for n=1:length(u.state)
                # calculate as sum of the momenta of the occupied
                get_momenta(du)[j,i] = get_momenta(du)[j,i] -sim.calculator.adiabatic_derivative[j,i][n, n]*u.state[n]
            end
        end
    end
=======
function acceleration!(dv, v, r, sim::Simulation{<:IESH}, t, state)
    #energ=zeros(1)
    for i=1:length(sim.atoms)
        for j=1:sim.DoFs
            #energ[1]=0
            #println(" ")
            dv[j,i] = 0.0
            for n=1:length(state)
                # calculate as sum of the momenta of the occupied
                dv[j,i] = dv[j,i] -sim.calculator.adiabatic_derivative[j,i][n, n]*state[n]
                #energ[1] = energ[1] + sim.calculator.eigenvalues[n]*u.state[n]
                #println(sim.calculator.eigenvalues[n]*u.state[n])
                                       #-sim.calculator.adiabatic_derivative[j,i][u.state, u.state]
            end
            dv[j,i] /= sim.atoms.masses[i]
        end
    end
    #println(energ)
>>>>>>> c0f6d97171f30e24cabdbf5157e3fd907556a7a1
end

# Calculate the _time_-derivate  of the  density matrix
# Goes back to motion!
# This one should be the same as for FSSH
<<<<<<< HEAD
function set_density_matrix_derivative!(du::IESHPhasespace, u::IESHPhasespace, sim::Simulation{<:IESH})
    σ = get_density_matrix(u)
    velocity = get_positions(du)
    V = sim.method.density_propagator # this is the potential energy surface

    V .= diagm(sim.calculator.eigenvalues)# Creates a diagonal matrix from eigenvalues
    #println(typeof(V))
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
    
        # 8 allocations, but still on average slightly  faster
    #@time 
    get_density_matrix(du) .= -im*(V*σ - σ*V)
    #inlining?
    # 4 allocations, but marginally slower on average
    #@time begin
    #mul!(sim.method.tmp_matrix1, σ, V)
    #mul!(sim.method.tmp_matrix2, V, σ)
    #mul!(sim.method.tmp_matrix1, get_density_matrix(u), V)
    #mul!(sim.method.tmp_matrix2, V, get_density_matrix(u))
    #get_density_matrix(du) .= -im *(sim.method.tmp_matrix2 - sim.method.tmp_matrix1)
    # No allocation, but slower
    #get_density_matrix(du) .= -im .*(sim.method.tmp_matrix2 .- sim.method.tmp_matrix1)
    #end
 end
=======
# function set_density_matrix_derivative!(du::IESHPhasespace, u::IESHPhasespace, sim::Simulation{<:IESH})
#     σ = get_density_matrix(u)
#     velocity = get_positions(du)
#     V = sim.method.density_propagator # this is the potential energy surface

#     V .= diagm(sim.calculator.eigenvalues)# Creates a diagonal matrix from eigenvalues
#     #println(typeof(V))
#     # Calculation going on here from: Martens_JPhysChemA_123_1110_2019, eq. 6
#     # d is the nonadiabatic coupling matrix
#     #i ħ dσ/dt = iħ sum_l [(V_{m,l} - i ħ velocity d_{m,l})*σ_{l,n} - &
#     #                      σ_{m,l}*(V_{l,n} - iħ velocity d_{l,n})]
#     # vgl. also SubotnikBellonzi_AnnuRevPhyschem_67_387_2016, eq. 5 and 6
#     for i=1:length(sim.atoms)
#         for j=1:sim.DoFs
#             V .-= im*velocity[j,i].*sim.method.nonadiabatic_coupling[j,i]
#         end
#     end
    
#     # This version is presumably much slower than below
#     #@time get_density_matrix(du) .= -im*(V*σ - σ*V)
#     #inlining?
#     ust = length(u.state)
#     Y = zeros(Complex, ust, ust)
#     Z = zeros(Complex, ust, ust)
#     # Okay, the loops are much slower than writing it with *
#     # Unfortunatley, this doesn't seem to speed things up
#     #@time begin
#     mul!(Y, σ, V)
#     mul!(Z, V, σ)
#     get_density_matrix(du) .= -im *(Z - Y)
#     #end
#  end
>>>>>>> c0f6d97171f30e24cabdbf5157e3fd907556a7a1

# This sets when the condition will be true. In this case, always.
condition(u, t, integrator::DiffEqBase.DEIntegrator) = true

function affect!(integrator::DiffEqBase.DEIntegrator)
<<<<<<< HEAD
    #println("ping7")
    #@time 
=======
>>>>>>> c0f6d97171f30e24cabdbf5157e3fd907556a7a1
    update_hopping_probability!(integrator)
    
    # not necessary anymore
    #new_state = select_new_state(integrator.p.method.hopping_probability, integrator.u.state)
    
    #if new_state != 0
    if integrator.p.method.hopping_probability[1] !=0
        #tick()
        println("Hop! from ", Int(integrator.p.method.hopping_probability[2]), " to ", Int(integrator.p.method.hopping_probability[3]))
        # Set new state population
        new_state = copy(integrator.u.state)
        new_state[Int(integrator.p.method.hopping_probability[2])] = 0
        new_state[Int(integrator.p.method.hopping_probability[3])] = 1


        # needs probably a vector as input for new_state (i.e. the state distribution)
        if calculate_rescaling_constant!(integrator, new_state)
            execute_hop!(integrator, new_state)
        else
            println("Frustrated!")
        end
        #tock()
    end
end

function update_hopping_probability!(integrator::DiffEqBase.DEIntegrator)
    sim = integrator.p
<<<<<<< HEAD
    coupling = sim.method.nonadiabatic_coupling
=======
    coupling = sim.calculator.nonadiabatic_coupling
>>>>>>> c0f6d97171f30e24cabdbf5157e3fd907556a7a1
    velocity = get_positions(get_du(integrator))
    s = integrator.u.state
    σ = get_density_matrix(integrator.u)
    dt = get_proposed_dt(integrator)
    
    sim.method.hopping_probability .= 0 # Set all entries to 0
    hop_mat = zeros(length(s),length(s))
    sumer=0
    sum_before = 0
    random_number = rand()
    first = true

    # Calculate matrix with hopping probs
    for l = 1:length(s)
        # If occupied
        if(integrator.u.state[l]==1)
            for m = 1:length(s)
                # if unoccupied
                if(integrator.u.state[m]==0)
                    for i=1:length(sim.atoms)
                        for j=1:sim.DoFs
                            hop_mat[l,m] += 2*velocity[j,i]*real(σ[m,l]/σ[l,l])*coupling[j,i][l,m] * dt
                            #println(velocity[j,i], " ",real(σ[m,l]/σ[l,l]), " ",coupling[j,i][l,m] , " ", dt)
                        end
                    end
                end    
                clamp(hop_mat[l,m], 0, 1)
                # Calculate the hopping probability. Hopping occures for
                # the transition that's first above the random number.
                # See: Tully_JChemPhys_93_1061_1990
                sumer = sumer+abs(hop_mat[l,m]) # cumulative sum.
                if (random_number > sumer &&first)
                    sum_before = sumer
                elseif (random_number < sumer && random_number > sum_before && first)
                    sim.method.hopping_probability[1] = sumer
                    sim.method.hopping_probability[2] = l
                    sim.method.hopping_probability[3] = m
                    first = false
                elseif (sumer > 1 && first)
                    println("Warning: Sum of hopping probability above 1!")
                    println("Sum: ", sumer, " Individ. hopping probability: ", hop_mat[l,m])
                    println("l = ", l, " m = ", m)
                    println(first)
                end
            end
        end
        #if (first == false) break
    end
    
 
end

# Calculate if hop frustrated
# This includes a momentum rescaling (IESH paper of Shenvi et al. expresses it in terms of velocity)
# I believe (also according to Subotnic&Miao JCP 150 2019) that this is the same.
# In any case, it's used to conserve energy. Reini remarked that we might not want it 
# eventually, but I'm leaving it for now.
# It should be related to: HammesSchifferTully_JChemPhys_101_4657_1994
function calculate_rescaling_constant!(integrator::DiffEqBase.DEIntegrator, new_state)::Bool
    sim = integrator.p
    old_state = integrator.u.state
    state_diff = integrator.p.method.hopping_probability
    velocity = get_positions(get_du(integrator))

    # Loop over and sum over optential energies, according to Eq. 12 Shenvi, Roy,  Tully,J. Chem. Phys. 130, 174107 (2009)
    c = 0
    for i=1:length(old_state)
        c = c + calculate_potential_energy_change(sim.calculator.eigenvalues, 
                                                  new_state[i], old_state[i],i)
    end
    a = zeros(length(sim.atoms))
    b = zero(a)
    # view: treats data structure from array as another array
    #': conjucated transposition (adjoint)
    @views for i in range(sim.atoms)
<<<<<<< HEAD
        coupling = [sim.method.nonadiabatic_coupling[j,i][Int(state_diff[3]), Int(state_diff[2])] for j=1:sim.DoFs]
=======
        coupling = [sim.calculator.nonadiabatic_coupling[j,i][Int(state_diff[3]), Int(state_diff[2])] for j=1:sim.DoFs]
>>>>>>> c0f6d97171f30e24cabdbf5157e3fd907556a7a1
        #coupling = [sim.method.nonadiabatic_coupling[j,i][Int(state_diff[2]), Int(state_diff[3])] for j=1:sim.DoFs]
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

# This does to update_electronics above.
function calculate_potential_energy_change(eigenvalues::Vector, new_state::Integer, 
    current_state::Integer, counter::Integer)
    new_state*eigenvalues[counter] - current_state*eigenvalues[counter]
end

function execute_hop!(integrator::DiffEqBase.DEIntegrator, new_state::Vector)
    state_diff = integrator.p.method.hopping_probability
    # For momentum rescaling, see eq. 7 and 8 SubotnikBellonzi_AnnuRevPhyschem_67_387_2016
    for i in range(integrator.p.atoms)
<<<<<<< HEAD
        coupling = [integrator.p.method.nonadiabatic_coupling[j,i][Int(state_diff[3]), 
=======
        coupling = [integrator.p.calculator.nonadiabatic_coupling[j,i][Int(state_diff[3]), 
>>>>>>> c0f6d97171f30e24cabdbf5157e3fd907556a7a1
                   Int(state_diff[2])] for j=1:integrator.p.DoFs]
        
        # No sum over states necessary, because momenta just rescaled & 
        # not calculated from scratch
<<<<<<< HEAD
        get_momenta(integrator.u)[:,i] .-= integrator.p.method.momentum_rescale[i] .* coupling
=======
        get_velocities(integrator.u)[:,i] .-= integrator.p.method.momentum_rescale[i] .* coupling ./ integrator.p.atoms.masses[i]
>>>>>>> c0f6d97171f30e24cabdbf5157e3fd907556a7a1

    end
    integrator.u.state = new_state
end