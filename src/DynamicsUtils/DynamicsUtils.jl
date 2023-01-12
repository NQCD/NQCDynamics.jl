
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

function classical_hamiltonian(sim, u) end

function classical_kinetic_energy(sim::AbstractSimulation, u)
    classical_kinetic_energy(sim,  DynamicsUtils.get_velocities(u))
end

function classical_kinetic_energy(sim::Simulation, v::AbstractMatrix)
    kinetic = zero(eltype(v))
    for i in axes(v, 2)
        for j in axes(v, 1)
            kinetic += masses(sim, i) * v[j,i]^2
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
    sum(V) + RingPolymers.get_spring_energy(sim.beads, masses(sim), r)
end

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
