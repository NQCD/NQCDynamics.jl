
"""
    DynamicsUtils

Utilities for dynamics simulations.
Includes:
* Basic dynamics variables functions
* Density matrix dynamics functions
* Standard callbacks to use during dynamics
* Plotting recipes for outputs
"""
module DynamicsUtils

using NQCDynamics:
    AbstractSimulation,
    Simulation,
    RingPolymerSimulation,
    Calculators,
    masses
using NQCDistributions: NonEqState
using Distributions: Bernoulli
using DelimitedFiles

"""
    divide_by_mass!(dv, masses)

Divide the contents of `dv` by the `masses`.
Assumes `dv` is an array of size (dofs, atoms) or (dofs, atoms, beads).
`masses` is the vector of masses for each atom that matches length with the second dimension.
"""
divide_by_mass!(dv, masses) = dv ./= masses'

multiply_by_mass!(dv, masses) = dv .*= masses'

"""
    velocity!(dr, v, r, sim, t)

Write the velocity `v` into `dr`.
Has extra arguments to work with `Dynamical(O/S)DEProblem`s.
"""
velocity!(dr, v, r, sim, t) = dr .= v

function acceleration! end

"""
    apply_interbead_coupling!(du::DynamicalVariables, u::DynamicalVariables,
                              sim::RingPolymerSimulation)
    
Applies the force that arises from the harmonic springs between adjacent beads.

Only applies the force for atoms labelled as quantum within the `RingPolymerParameters`.
"""
function apply_interbead_coupling!(dr::AbstractArray{T,3}, r::AbstractArray{T,3}, sim::RingPolymerSimulation) where {T}
    for i in axes(dr, 3)
        iplus = mod1(i+1, sim.beads.n_beads)
        iminus = mod1(i-1, sim.beads.n_beads)
        for j in sim.beads.quantum_atoms
            for k in axes(dr, 1)
                dr[k,j,i] -= sim.beads.ω_n² * (2r[k,j,i] - r[k,j,iplus] - r[k,j,iminus])
            end
        end
    end
    return nothing
end

function classical_hamiltonian(sim::Simulation, u)
    kinetic = classical_kinetic_energy(sim, u)
    potential = classical_potential_energy(sim, u)
    return kinetic + potential
end

function classical_hamiltonian(sim::RingPolymerSimulation, u)
    spring = classical_spring_energy(sim, u)
    kinetic = classical_kinetic_energy(sim, u)
    potential = classical_potential_energy(sim, u)
    return spring + kinetic + potential
end

function centroid_classical_hamiltonian(sim::RingPolymerSimulation, u)
    kinetic = centroid_classical_kinetic_energy(sim, u)
    potential = centroid_classical_potential_energy(sim, u)
    return kinetic + potential
end

function centroid_classical_kinetic_energy(sim::RingPolymerSimulation, u)
    v = DynamicsUtils.get_velocities(u)
    centroid_v = get_centroid(v)
    kinetic = DynamicsUtils.classical_kinetic_energy(masses(sim), centroid_v)
    return kinetic
end

function centroid_classical_potential_energy end

function classical_spring_energy(sim::RingPolymerSimulation, u)
    return classical_spring_energy(sim, get_positions(u))
end

function classical_spring_energy(sim::RingPolymerSimulation, r::AbstractArray{T,3}) where {T}
    return RingPolymers.get_spring_energy(sim.beads, sim.atoms.masses, r)
end

function classical_kinetic_energy(sim::AbstractSimulation, u)
    classical_kinetic_energy(sim,  DynamicsUtils.get_velocities(u))
end

function classical_kinetic_energy(sim::AbstractSimulation, v::AbstractMatrix)
    return classical_kinetic_energy(masses(sim), v)
end

function classical_kinetic_energy(mass::AbstractVector, v::AbstractMatrix)
    kinetic = zero(eltype(v))
    for i in axes(v, 2)
        for j in axes(v, 1)
            kinetic += mass[i] * v[j,i]^2
        end
    end
    return kinetic / 2
end

function classical_kinetic_energy(sim::RingPolymerSimulation, v::AbstractArray{T,3}) where {T}
    kinetic = zero(eltype(v))
    for k in axes(v, 3)
        for i in axes(v, 2)
            for j in axes(v, 1)
                kinetic += masses(sim, i) * v[j,i,k]^2
            end
        end
    end
    return kinetic / 2
end

function classical_potential_energy(sim::AbstractSimulation, u)
    classical_potential_energy(sim,  DynamicsUtils.get_positions(u))
end

function classical_potential_energy(sim::Simulation, r::AbstractMatrix)
    Calculators.get_potential(sim.calculator, r)
end

function classical_potential_energy(sim::RingPolymerSimulation, r::AbstractArray{T,3}) where {T}
    V = Calculators.get_potential(sim.calculator, r)
    return sum(V)
end

function get_hopping_eigenvalues end
function get_hopping_nonadiabatic_coupling end
function get_hopping_velocity end

include("dynamics_variables.jl")
include("callbacks.jl")
export CellBoundaryCallback
export TerminatingCallback

include("density_matrix_dynamics.jl")
include("wavefunction_dynamics.jl")
include("plot.jl")

function set_unoccupied_states!(unoccupied::AbstractVector, occupied::AbstractVector)
    nstates = length(unoccupied) + length(occupied) 
    index = 1
    for i in 1:nstates
        if !(i in occupied)
            unoccupied[index] = i
            index += 1
        end
    end
end

fermi(ϵ, μ, β) = 1 / (1 + exp(β*(ϵ-μ)))

function sample_fermi_dirac_distribution(energies, nelectrons, available_states, β)
    nstates = length(available_states)
    state = collect(Iterators.take(available_states, nelectrons)) # populate your states according to the number of electrons you have
    for _ in 1:(nstates * nelectrons) # iterate many times where you check to see if you should make a state change
        current_index = rand(eachindex(state)) # makes rand() return a random index of state array instead of a random element
        i = state[current_index] # Pick random occupied state
        j = rand(setdiff(available_states, state)) # Pick random unoccupied state
        prob = exp(-β * (energies[j] - energies[i])) # Calculate Boltzmann factor - replace this with bernouili of non-eq dist
        if prob > rand()
            state[current_index] = j # Set unoccupied state to occupied
        end
    end
    sort!(state)
    return state
end

# ----------------------------------------- NEW FUNCTIONS ---------------------------------------- #
# Insert here the wrapper for Henry's DiscretizationScript - may want to move this to NQCDistributions instead?

include("noneqdist_discretization.jl")

"""
    NonEqDist(dist_filename::String, DOS_filename::String)

Non-equilibrium distribution and DoS are read in from file and mapped to the energy grid generated for the system.
The `NonEqState` electronic distribution is output.
"""
function NonEqDist(dist_filename::String, DOS_filename::String)

    @info "reading in from external files"
    dis_from_file = readdlm(dist_filename) # Read distribution from file
    energy_grid = dis_from_file[:,1] # Seperate into energy grid 
    distribution = dis_from_file[:,2] # and total distribution
    dis_spl = LinearInterpolation(distribution,energy_grid) # Spline the distribution to project onto new energy grid
    DOS_spl = generate_DOS(DOS_filename,85) # Build DOS spline from file - n = 85 corresponds to the interpolation factor

    @info "create `NonEqState` electronic state"
    return NonEqState(dis_spl, DOS_spl)
end

"""
    NonEqDist(sample_distribution::Vector{Float64})

Takes existing non-equilibrium distribution and energy grid in array format and converts to spline functions.
(DOS is not included in this version, `dis_spl` is mapped to both fields temporarily!)
"""
function NonEqDist(sample_distribution::Vector{Float64}, sample_energy_grid::Vector{Float64})

    dis_spl = LinearInterpolation(sample_distribution, sample_energy_grid; extrapolate=true) # Spline the distribution to project onto new energy grid
    # DOS_spl = generate_DOS(DOS_filename,85) # Build DOS spline from file - n = 85 corresponds to the interpolation factor

    @info "create `NonEqState` electronic state"
    return NonEqState(dis_spl, dis_spl)
end


"""
    sample_noneq_distribution(distribution, nelectrons, available_states)

New sampling function which generates a state object from a pre-existing non-equilibrium distribution.
Distribution must have been mapped to the number of available states.
"""
function sample_noneq_distribution(energies, nelectrons, available_states, dis_spline)

    # DOS = DOS_spline(energies) # Recast DOS and distribution onto new grid
    dis = dis_spline(energies)

    # confining values in dis to have upper bounds of 0.0 and 1.0
    dis[dis .< 0.0] .= 0.0
    dis[dis .> 1.0] .= 1.0

    # add check here to ensure distribution is same dimension as available_states
    nstates = length(available_states)
    state = collect(Iterators.take(available_states, nelectrons)) # populate your states according to the number of electrons you have
    for _ in 1:(nstates * nelectrons) # iterate many times where you check to see if you should make a state change
        current_index = rand(eachindex(state)) # makes rand() return a random index of state array instead of a random element
        i = state[current_index] # Pick random occupied state
        prob_RemainOccupied = Bernoulli(dis[i]) # Bernoulli of non-eq distribution at selected occupied state
        if !rand(prob_RemainOccupied) # check if occupied state should remain occupied
            j = rand(setdiff(available_states, state)) # Pick random unoccupied state
            prob_BecomeOccupied = Bernoulli(dis[j]) # Bernouili of non-eq distribution at selected unoccupied state
            if rand(prob_BecomeOccupied) # will return True if a random number is less than the probability given by Bernoulli
                state[current_index] = j # Set unoccupied state to occupied
            end
        end
    end
    sort!(state)
    return state
end

### TEMP ###
function DEBUG_sample_noneq_distribution(energies, nelectrons, available_states, dis_spline)

    # DOS = DOS_spline(energies) # Recast DOS and distribution onto new grid
    dis = dis_spline(energies)

    # counters for returning proportion successful
    count_UnoccupySource = 0
    count_OccupyDestination = 0

    # add check here to ensure distribution is same dimension as available_states
    nstates = length(available_states)
    total_iterations = nstates * nelectrons
    state = collect(Iterators.take(available_states, nelectrons)) # populate your states according to the number of electrons you have
    for _ in 1:(total_iterations) # iterate many times where you check to see if you should make a state change
        current_index = rand(eachindex(state)) # makes rand() return a random index of state array instead of a random element
        i = state[current_index] # Pick random occupied state
        prob_RemainOccupied = Bernoulli(dis[i]) # Bernoulli of non-eq distribution at selected occupied state
        if !rand(prob_RemainOccupied) # check if occupied state should remain occupied
            count_UnoccupySource += 1
            j = rand(setdiff(available_states, state)) # Pick random unoccupied state
            prob_BecomeOccupied = Bernoulli(dis[j]) # Bernouili of non-eq distribution at selected unoccupied state
            if rand(prob_BecomeOccupied) # will return True if a random number is less than the probability given by Bernoulli
                count_OccupyDestination += 1
                state[current_index] = j # Set unoccupied state to occupied
            end
        end
    end
    sort!(state)
    return state, (total_iterations, count_UnoccupySource, count_OccupyDestination)
end
### TEMP ###
"""
    simple_non_eq(ϵ, μ, β, Δf_β1, Δf_β2)

Outputs a simple approximation to a non-equilibrium distribution, generated from finite temperature Fermi-Dirac distributions, that can be used for test purposes.
"""
function simple_non_eq(ϵ, μ, β, Δf_β1, Δf_β2)
    β != Inf || throw(error("Non-equilibrium distribution needs to be a finite temperature, β cannot be Inf."))
    Δf = fermi.(ϵ, μ, Δf_β1) - fermi.(ϵ, μ, Δf_β2) #  create perterbation from 2 finite temperature Fermi-Dirac distributions
    F = fermi.(ϵ, μ, β) + Δf  # add perterbation to a finite temperature Fermi-Dirac distribution
    return F
end

# ------------------------------------------------------------------------------------------------ #

get_available_states(::Colon, nstates::Integer) = 1:nstates
function get_available_states(available_states::AbstractVector, nstates::Integer)
    maximum(available_states) > nstates && throw(DomainError(available, "There are only $nstates in the system."))
    return available_states
end


end # module
