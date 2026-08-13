
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
    masses
using NQCCalculators

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
function velocity!(dr, v, r, sim, t)
    @inbounds for i in eachindex(dr)
        dr[i] = v[i]
    end
    return nothing
end

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
    NQCCalculators.update_potential!(sim.cache, r)
    return NQCCalculators.get_potential(sim.cache, r)[1]
end

function classical_potential_energy(sim::RingPolymerSimulation, r::AbstractArray{T,3}) where {T}
    NQCCalculators.update_potential!(sim.cache, r)
    V = NQCCalculators.get_potential(sim.cache, r)
    return sum(V)[1]
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


function sample_fermi_dirac_distribution(energies, nelectrons, available_states, β, μ) #Uses ratio of Fermi probabilities
    nstates = length(available_states)
    state = collect(Iterators.take(available_states, nelectrons))
    for _ in 1:(nstates * nelectrons)
        current_index = rand(eachindex(state))
        i = state[current_index] # Pick random occupied state
        j = rand(setdiff(available_states, state)) # Pick random unoccupied state
        f1 = 1/(exp((energies[i]-μ)β)+1)
        f2 = 1/(exp((energies[j]-μ)β)+1)
        prob = f2/f1 # Calculate Boltzmann factor
        if prob > rand()
            state[current_index] = j # Set unoccupied state to occupied
        end
    end
    sort!(state)
    return state
end

function sample_fermi_dirac_distribution(energies, nelectrons, available_states, β) # Uses Boltzmann Factor
    nstates = length(available_states)
    state = collect(Iterators.take(available_states, nelectrons))
    for _ in 1:(nstates * nelectrons)
        current_index = rand(eachindex(state))
        i = state[current_index] # Pick random occupied state
        j = rand(setdiff(available_states, state)) # Pick random unoccupied state
        prob = exp(-β * (energies[j] - energies[i])) # Calculate Boltzmann factor
        if prob > rand()
            state[current_index] = j # Set unoccupied state to occupied
        end
    end
    sort!(state)
    return state
end

get_available_states(::Colon, nstates::Integer) = 1:nstates
function get_available_states(available_states::AbstractVector, nstates::Integer)
    maximum(available_states) > nstates && throw(DomainError(available, "There are only $nstates in the system."))
    return available_states
end


end # module
