using StatsBase: mean
using .Calculators: DiabaticCalculator, RingPolymerDiabaticCalculator

export wave_IESH
export SurfaceHoppingVariablesIESH
export get_wavefunction_matrix

"""This module controles how IESH is executed. For a description of IESH, see e.g.
Roy, Shenvi, Tully, J. Chem. Phys. 130, 174716 (2009) and 
Shenvi, Roy, Tully, J. Chem. Phys. 130, 174107 (2009).

The density matrix is set up in surface_hopping_variables.jl
"""

abstract type SurfaceHoppingIESH <: Method end

mutable struct wave_IESH{T} <: SurfaceHoppingIESH
    wavefunction_propagator::Matrix{Complex{T}}
    hopping_probability::Vector{T}
    new_state::Vector{Int}
    function wave_IESH{T}(states::Integer) where {T}

    
        wavefunction_propagator = zeros(states, states)
        hopping_probability = zeros(3)
        new_state = zeros(states) # this probably needs to be modified
        new{T}(wavefunction_propagator, hopping_probability, new_state)
    end
end

# Define for wave function propagation in IESH. It might be possible to unite this w/ the top
mutable struct SurfaceHoppingVariablesIESH{T,D,S}  <: DynamicalVariables{T}
    x::ArrayPartition{Complex{T}, Tuple{D,D,Matrix{Complex{T}}}}
    state::S
end

function SurfaceHoppingVariablesIESH(x::ArrayPartition{T}, state) where {T<:AbstractFloat}
    SurfaceHoppingVariablesIESH(ArrayPartition(x.x[1], x.x[2], Complex.(x.x[3])), state)
end

"""Set up matrix that stores wave vectors for propagation
See ShenviRoyTully_JChemPhys_130_174107_200
Note that different from the density matrix, wave_mat defines vectors for a state, 
where 1 indicated the state, but not whether occupied or not."""
function SurfaceHoppingVariablesIESH(v::AbstractArray, r::AbstractArray, n_states::Integer, state::Vector{Int})
    wave_mat = zeros(Complex{eltype(r)}, n_states, n_states)
    #for i=1:n_states/2
    for i=1:n_states
        wave_mat[Int(i), Int(i)] = 1
    end
    SurfaceHoppingVariablesIESH(ArrayPartition(v, r, wave_mat), state)
end

get_wavefunction_matrix(u::SurfaceHoppingVariablesIESH) = u.x.x[3]

"""motion! is given to ODE Problem and propagated by the timestep there"""
function motion!(du, u, sim::AbstractSimulation{<:SurfaceHoppingIESH}, t)
    #println("ping2")
    dr = get_positions(du)
    dv = get_velocities(du)
    dσ = get_wavefunction_matrix(du)


    r = get_positions(u)
    v = get_velocities(u)
    σ = get_wavefunction_matrix(u)
    println(dσ[1,1], σ[1,1])

    # presumably comes from DifferentialEquations-Julia module
    velocity!(dr, v, r, sim, t)
    # src/Calculators/Calculators.jl
    # Gets energy, forces, eigenvalues and nonadiabatic couplings.
    Calculators.update_electronics!(sim.calculator, r) # uses only nuclear DOF
    acceleration!(dv, v, r, sim, t; state=u.state) # nuclear DOF
    set_wavefunction_derivative!(dσ, v, σ, sim, u)
end


"""This is part of the adiabatic propagation of the nuclei.
   See Eq. 12 of Shenvi, Tully JCP 2009 paper."""
function acceleration!(dv, v, r, sim::Simulation{<:wave_IESH}, t; state=1)
    #println("ping3")
    # This seems to look all right
    # Goes over direction 2
    dv .= 0.0
    for i in axes(dv, 2)
        for j in axes(dv, 1)        
            for k in 1:length(state)
                # Include only occupied state.
                # It's kind of double with *state[k].
                if state[k] == 1
                    # Calculate as the sum of the momenta of occupied states
                    dv[j,i] = dv[j,i] - sim.calculator.adiabatic_derivative[j,i][k, k]*
                              state[k] / sim.atoms.masses[i]
                end
            end
        end
    end
    return nothing
end

"""Propagation of electronic wave function happens according to Eq. (14) 
   in the Shenvi, Tully paper (JCP 2009)
   The extended formula is taken from HammesSchifferTully_JChemPhys_101_4657_1994, Eq. (16):
   iħ d ψ_{k}/dt = ∑_j ψ_{j}(V_{kj} - i v d_{kj})
   Where v is the velocity. See also Tully_JChemPhys_93_1061_1990, Eq. (7). 
   According to point (3) of the algorithm of Shenvi, Tully, 2009, only the occupied
   wave functions are integrated. j, I believe, should however run other the unoccupied
   nonetheless, because still, even if not propagated, all wavefunctions should contribute
   to dσ (also, otherwise, no hopping probability)
   """
function set_wavefunction_derivative!(dσ, v, σ, sim::Simulation{<:SurfaceHoppingIESH}, u)
    #println("ping4")
    V = sim.method.wavefunction_propagator
    s = u.state

    # Get the eigenvalues
    # Electronic H.
    V .= diagm(sim.calculator.eigenvalues)
    ## Try the diabatic representation of V
    ##V1 = zeros(length(s),length(s))
    ##V1 .= diagm(sim.calculator.eigenvalues)
    ##V .= sim.calculator.eigenvectors * V1 * sim.calculator.eigenvectors'
    n_states = sim.calculator.model.n_states
    dσ .= 0.0
    # Inner loop over j may not be necessary, since V_{j,i} should be zero.
    for k = 1:n_states
        # is occupied?
        if (s[k] == 1)
            for j = 1:n_states
                # occupied?
                #if (s[j] == 1)
                    for I in eachindex(v)
                        dσ[k,:] .= dσ[k,:] + σ[j,:]*(V[j,k] - im*v[I]*
                                   sim.calculator.nonadiabatic_coupling[I][j,k])
                                   bbb = σ[j,:]*(V[j,k] - im*v[I]*
                                   sim.calculator.nonadiabatic_coupling[I][j,k])
                        #println(k, " ", j," ",σ[j,:]," ",sim.calculator.nonadiabatic_coupling[I][j,k], " ", V[j,k], " ", real(bbb))
                    end
                #end
            end
        end
    end
end

function check_hop!(u, t, integrator)::Bool
    #println("ping5")
    evaluate_hopping_probability!(
        integrator.p,
        u,
        get_proposed_dt(integrator))

    integrator.p.method.new_state = select_new_state(integrator.p, u)
    return integrator.p.method.new_state != u.state
end

"""Hopping probability according to Eq.s (21), (17), (18) in Shenvi/Roy/Tully JCP 2009 paper.
   The density matrix is used there, so I am setting up the density matrix here, too"""
function evaluate_hopping_probability!(sim::Simulation{<:wave_IESH}, u, dt)
    #println("ping6")
    v = get_velocities(u)
    Ψ = get_wavefunction_matrix(u)
    s = u.state
    d = sim.calculator.nonadiabatic_coupling
    n_states = sim.calculator.model.n_states

    hop_mat = zeros(n_states, n_states)
    sim.calculator.tmp_mat_complex1 .= 0
    sumer = 0
    first = true
    random_number = rand()
    sum_before = 0.0
    lc = 0
    mc = 0

    sim.method.hopping_probability .= 0 # Set all entries to 0
    # Assemble density matrix (multiplied by occuplation)
    for i = 1:n_states
        sim.calculator.tmp_mat_complex1 .+= (Ψ[i,:] * Ψ[i,:]')*s[i]
    end
    #println(real(sim.calculator.tmp_mat_complex1))

    for l = 1:n_states
        # Is occupied?
        if(s[l] == 1)
            for m = 1:n_states
                # Is unoccupied?
                if (s[m] == 0)
                    for I in eachindex(v)
                        #sim.method.hopping_probability[m] += 2v[I]*real(σ[m,s]/σ[s,s])*d[I][s,m] * dt
                        hop_mat[l,m] = 2*v[I]*real(sim.calculator.tmp_mat_complex1[m,l]/
                                                   sim.calculator.tmp_mat_complex1[l,l])*d[I][l,m] * dt
                        #println(v[I], " ", real(sim.calculator.tmp_mat_complex1[m,l]), " ",
                        #real(sim.calculator.tmp_mat_complex1[l,l]), " ", d[I][l,m], " ", dt)
                        #println(real(sim.calculator.tmp_mat_complex1[l,m]))
                    end
                end # end if 
                clamp(hop_mat[l,m], 0, 1)
                # Calculate the hopping probability. Hopping occures for
                # the transition that's first above the random number.
                # See: Tully_JChemPhys_93_1061_1990
                sumer = sumer + abs(hop_mat[l,m]) # cumulative sum.
                #println(sumer)
                # If sum of hopping probabilities is larger than random number,
                # hopping can occur
                if (random_number > sumer && first)
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
    end
    return nothing
end

"""
Set up new states for hopping
"""
function select_new_state(sim::AbstractSimulation{<:wave_IESH}, u)
    #println("ping7")
    if sim.method.hopping_probability[1] !=0
        println("Hop! from ", Int(sim.method.hopping_probability[2]), " to ", Int(sim.method.hopping_probability[3]))
        # Set new state population
        new_state = copy(u.state)
        new_state[Int(sim.method.hopping_probability[2])] = 0
        new_state[Int(sim.method.hopping_probability[3])] = 1
        return new_state
    end
    return u.state
end

function execute_hop!(integrator)
    #println("ping8")
    rescale_velocity!(integrator.p, integrator.u) && (integrator.u.state = integrator.p.method.new_state)
    return nothing
end

function rescale_velocity!(sim::AbstractSimulation{<:wave_IESH}, u)::Bool
    #println("ping9")
    old_state = u.state
    new_state = sim.method.new_state
    state_diff = sim.method.hopping_probability
    velocity = get_velocities(u)
    
    # Loop over and sum over potential energies, according to Eq. 12 Shenvi, Roy,  Tully,J. Chem. Phys. 130, 174107 (2009)
    c = 0
    for i=1:length(old_state)
        c = c + calculate_potential_energy_change(sim.calculator, 
                                                  new_state[i], old_state[i], i)
    end
    a, b = evaluate_a_and_b(sim, velocity, state_diff)
    discriminant = b.^2 .- 2a.*c
    
    any(discriminant .< 0) && println("Frustrated!")
    any(discriminant .< 0) && return false

    root = sqrt.(discriminant)
    velocity_rescale = min.(abs.((b .+ root) ./ a), abs.((b .- root) ./ a))
    perform_rescaling!(sim, velocity, velocity_rescale, state_diff)

    return true
end

"""
Evaluate nonadiabatic coupling after hop
"""
function evaluate_a_and_b(sim::AbstractSimulation{<:wave_IESH}, velocity::AbstractArray, state_diff)
    #println("ping10")
    a = zeros(length(sim.atoms))
    b = zero(a)
    @views for i in range(sim.atoms)
        coupling = [sim.calculator.nonadiabatic_coupling[j,i][Int(state_diff[3]), 
                    Int(state_diff[2])] for j=1:sim.DoFs]
        a[i] = coupling'coupling / sim.atoms.masses[i]
        b[i] = velocity[:,i]'coupling
    end
    return (a, b)
end


"""
Performs momentum rescaling, see eq. 7 and 8 SubotnikBellonzi_AnnuRevPhyschem_67_387_2016
"""
function perform_rescaling!(sim::Simulation{<:wave_IESH}, velocity, velocity_rescale,
                            state_diff)
    #println("ping11")
    for i in range(sim.atoms)
        coupling = [sim.calculator.nonadiabatic_coupling[j,i][Int(state_diff[3]), 
                    Int(state_diff[2])] for j=1:sim.DoFs]
        velocity[:,i] .-= velocity_rescale[i] .* coupling ./ sim.atoms.masses[i]
    end
    return nothing
end

# Goes to rescale_velocity!
function calculate_potential_energy_change(calc::DiabaticCalculator, new_state::Integer, 
                                           current_state::Integer, counter::Integer)
    #println("ping12")
    #DeltaE = calc.eigenvalues[counter]*new_state - calc.eigenvalues[counter]*current_state
    return calc.eigenvalues[counter]*new_state - calc.eigenvalues[counter]*current_state
end


const HoppingCallback = DiscreteCallback(check_hop!, execute_hop!; save_positions=(false, false))
get_callbacks(::AbstractSimulation{<:SurfaceHoppingIESH}) = HoppingCallback

"""
This function should set the field `sim.method.hopping_probability`.
"""
function evaluate_hopping_probability!(::AbstractSimulation{<:SurfaceHoppingIESH}, u, dt) end

"""
This function should return the desired state determined by the probability.
Should return the original state if no hop is desired.
"""
function select_new_state(::AbstractSimulation{<:SurfaceHoppingIESH}, u) end

"""
This function should modify the velocity and return a `Bool` that determines
whether the state change should take place.

This only needs to be implemented if the velocity should be modified during a hop.
"""
rescale_velocity!(::AbstractSimulation{<:SurfaceHoppingIESH}, u) = true

create_problem(u0, tspan, sim::AbstractSimulation{<:SurfaceHoppingIESH}) = 
               ODEProblem(motion!, u0, tspan, sim)

