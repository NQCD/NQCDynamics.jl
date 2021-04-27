using StatsBase: mean
using .Calculators: DiabaticCalculator, RingPolymerDiabaticCalculator

export IESH

"""This module controles how IESH is executed. For a description of IESH, see e.g.
Roy, Shenvi, Tully, J. Chem. Phys. 130, 174716 (2009) and 
Shenvi, Roy, Tully, J. Chem. Phys. 130, 174107 (2009).
"""

mutable struct IESH{T} <: SurfaceHopping
    density_propagator::Matrix{Complex{T}}
    hopping_probability::Vector{T}
    new_state::Vector{Int}
    function IESH{T}(states::Integer) where {T}
        density_propagator = zeros(states, states)
        hopping_probability = zeros(3)
        new_state = zeros(states) # this probably needs to be modified
        new{T}(density_propagator, hopping_probability, new_state)
    end
end

# See Eq. 12 of Shenvi, Tully paper.
function acceleration!(dv, v, r, sim::Simulation{<:IESH}, t; state=1)
    #println("ping6")
    # Goes over direction 2
    for i in axes(dv, 2)
        for j in axes(dv, 1)
            #sim.calculator.adiabatic_derivative[j,i] .= 0.0
            for k in 1:length(state)
                # Calculate as the sum of the momenta of occupied states
                dv[j,i] = sim.calculator.adiabatic_derivative[j,i][k, k]-
                          sim.calculator.adiabatic_derivative[j,i][k, k]*
                          state[k] / sim.atoms.masses[i]
            end
        end
    end
    return nothing
end

function evaluate_hopping_probability!(sim::Simulation{<:IESH}, u, dt)
    #println("ping7")
    v = get_velocities(u)
    σ = get_density_matrix(u)
    s = u.state
    d = sim.calculator.nonadiabatic_coupling
    n_st = sim.calculator.model.n_states

    hop_mat = zeros(n_st, n_st)
    sumer = 0
    first = true
    random_number = rand()
    sum_before = sumer

    sim.method.hopping_probability .= 0 # Set all entries to 0
    for l = 1:n_st
        # Is occupied?
        if(s[l] == 1)
            for m = 1:n_st
                # Is unoccupied?
                if (s[m] == state)
                    for I in eachindex(v)
                        #sim.method.hopping_probability[m] += 2v[I]*real(σ[m,s]/σ[s,s])*d[I][s,m] * dt
                        hop_mat[l,m] = 2*v[I]*real(σ[m,l]/σ[l,l])*d[I][l,m] * dt
                    end
                end
                clamp(hop_mat[l,m], 0, 1)
                # Calculate the hopping probability. Hopping occures for
                # the transition that's first above the random number.
                # See: Tully_JChemPhys_93_1061_1990
                sumer = sumer + abs(hop_mat[l,m]) # cumulative sum.
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
    #clamp!(sim.method.hopping_probability, 0, 1) # Restrict probabilities between 0 and 1
    #cumsum!(sim.method.hopping_probability, sim.method.hopping_probability)
    return nothing
end

"""
Set up new states for hopping
"""
function select_new_state(sim::AbstractSimulation{<:IESH}, u)
    #println("ping8")
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

function rescale_velocity!(sim::AbstractSimulation{<:IESH}, u)::Bool
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

    any(discriminant .< 0) && println("Frustrated!") && return false

    root = sqrt.(discriminant)
    velocity_rescale = min.(abs.((b .+ root) ./ a), abs.((b .- root) ./ a))
    perform_rescaling!(sim, velocity, velocity_rescale, state_diff)

    return true
end

"""
Evaluate nonadiabatic coupling after hop
"""
function evaluate_a_and_b(sim::AbstractSimulation{<:IESH}, velocity::AbstractArray, state_diff)
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
function perform_rescaling!(sim::Simulation{<:IESH}, velocity, velocity_rescale,
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

