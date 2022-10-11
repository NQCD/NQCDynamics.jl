
"""
    DynamicsOutputs

Defines a set of functions that can be used to calculate outputs for dynamics simulations.
"""
module DynamicsOutputs

using RingPolymerArrays: get_centroid
using UnitfulAtomic: austrip
using LinearAlgebra: norm
using ComponentArrays: ComponentVector

using NQCDynamics:
    Estimators,
    DynamicsUtils

using NQCModels: NQCModels

using .DynamicsUtils:
    get_positions,
    get_velocities,
    get_quantum_subsystem

using ..InitialConditions: QuantisedDiatomic

"""
    OutputVelocity(sol, i)

Output the velocity at each timestep during the trajectory.
"""
OutputVelocity(sol, i) = [copy(get_velocities(u)) for u in sol.u]
export OutputVelocity

"""
    OutputCentroidVelocity(sol, i)

Output the velocity of the ring polymer centroid at each timestep during the trajectory. 
"""
OutputCentroidVelocity(sol, i) = [get_centroid(get_velocities(u)) for u in sol.u]
export OutputCentroidVelocity

"""
    OutputPosition(sol, i)

Output the position at each timestep during the trajectory.
"""
OutputPosition(sol, i) = [copy(get_positions(u)) for u in sol.u]
export OutputPosition

"""
    OutputCentroidPosition(sol, i)

Output the position of the ring polymer centroid at each timestep during the trajectory.
"""
OutputCentroidPosition(sol, i) = [get_centroid(get_position(u)) for u in sol.u]
export OutputCentroidPosition

"""
    OutputPotentialEnergy(sol, i)

Output the adiabatic potential energy at each timestep during the trajectory.
"""
OutputPotentialEnergy(sol, i) = DynamicsUtils.classical_potential_energy.(sol.prob.p, sol.u)
export OutputPotentialEnergy

"""
    OutputTotalEnergy(sol, i)

Evaluate the classical Hamiltonian at each timestep during the trajectory.
"""
OutputTotalEnergy(sol, i) = DynamicsUtils.classical_hamiltonian.(sol.prob.p, sol.u)
export OutputTotalEnergy

"""
    OutputKineticEnergy(sol, i)

Evaluate the classical kinetic energy at each timestep during the trajectory.
"""
OutputKineticEnergy(sol, i) = DynamicsUtils.classical_kinetic_energy.(sol.prob.p, sol.u)
export OutputKineticEnergy

"""
    OutputDynamicsVariables(sol, i)

Output all of the dynamics variables at each timestep during the trajectory.
"""
OutputDynamicsVariables(sol, i) = copy.(sol.u)
export OutputDynamicsVariables

"""
    OutputQuantumSubsystem(sol, i)

Output the quantum subsystem at each timestep during the trajectory.
Usually this will refer to a wavefunction or density matrix but will depend on the particular dynamics method.
"""
OutputQuantumSubsystem(sol, i) = [copy(get_quantum_subsystem(u)) for u in sol.u]
export OutputQuantumSubsystem

"""
    OutputMappingPosition(sol, i)

Output the position mapping variables at each timestep during the trajectory.
"""
OutputMappingPosition(sol, i) = [copy(DynamicsUtils.get_mapping_positions(u)) for u in sol.u]
export OutputMappingPosition

"""
    OutputMappingMomentum(sol, i)

Output the momentum mapping variable at each timestep during the trajectory.
"""
OutputMappingMomentum(sol, i) = [copy(DynamicsUtils.get_mapping_momenta(u)) for u in sol.u]
export OutputMappingMomentum

"""
    OutputDiscreteState(sol, i)

Output the discrete state variable at each timestep during the trajectory.
This is used for surface hopping simulations and returns the variable that determines the currently occupied adiabatic state.

Requires that the dynamics variable has a field `state`.

Use [`OutputDiabaticPopulation`](@ref) or [`OutputAdiabaticPopulation`](@ref) to get the population estimators.
"""
OutputDiscreteState(sol, i) = [copy(u.state) for u in sol.u]
export OutputDiscreteState

"""
    OutputDiabaticPopulation(sol, i)

Output the diabatic population at each timestep during the trajectory.
"""
OutputDiabaticPopulation(sol, i) = Estimators.diabatic_population.(sol.prob.p, sol.u)
export OutputDiabaticPopulation

"""
    OutputTotalDiabaticPopulation(sol, i)

Output the total diabatic population at eah timestep during the trajectory.
"""
OutputTotalDiabaticPopulation(sol, i) = sum.(Estimators.diabatic_population(sol.prob.p, sol.u))
export OutputTotalDiabaticPopulation

"""
    OutputAdiabaticPopulation(sol, i)

Output the adiabatic population at each timestep during the trajectory.
"""
OutputAdiabaticPopulation(sol, i) = Estimators.adiabatic_population.(sol.prob.p, sol.u)
export OutputAdiabaticPopulation

"""
    OutputTotalAdiabaticPopulation(sol, i)

Output the total adiabatic population at each timestep during the trajectory.
"""
OutputTotalAdiabaticPopulation(sol, i) = sum.(Estimators.adiabatic_population.(sol.prob.p, sol.u))
export OutputTotalAdiabaticPopulation

"""
Output the end point of each trajectory.
"""
OutputFinal(sol, i) = last(sol.u)
export OutputFinal

"""
Output a 1 if the molecule has dissociated, 0 otherwise.
"""
struct OutputDissociation{T}
    "The maximum distance at which the two atoms can be considered bonded."
    distance::T
    "The indices of the two atoms in the molecule of interest."
    atom_indices::Tuple{Int,Int}
    OutputDissociation(distance, atom_indices) = new{typeof(distance)}(austrip(distance), atom_indices)
end
export OutputDissociation

function (output::OutputDissociation)(sol, i)
    R = DynamicsUtils.get_positions(last(sol.u))
    dissociated = norm(R[:,output.atom_indices[1]] .- R[:,output.atom_indices[2]]) > output.distance
    return dissociated ? 1 : 0
end

"""
Output the vibrational and rotational quantum numbers of the final image.
"""
struct OutputQuantisedDiatomic{S,H,V}
    sim::S
    height::H
    normal_vector::V
end
OutputQuantisedDiatomic(sim; height=10, normal_vector=[0, 0, 1]) = OutputQuantisedDiatomic(sim, height, normal_vector)
export OutputQuantisedDiatomic

function (output::OutputQuantisedDiatomic)(sol, i)
    final = last(sol.u) 
    ν, J = QuantisedDiatomic.quantise_diatomic(output.sim,
        DynamicsUtils.get_velocities(final), DynamicsUtils.get_positions(final);
        height=output.height, normal_vector=output.normal_vector)
    return (ν, J)
end

"""
Output a `ComponentVector` with fields `reflection` and `transmission` containing
the probability of the outcome.
Each index in the arrays refers to the adiabatic state.
"""
struct OutputStateResolvedScattering1D{S}
    sim::S
    type::Symbol
end
function (output::OutputStateResolvedScattering1D)(sol, i)
    final = last(sol.u) # get final configuration from trajectory
    if output.type == :adiabatic
        populations = Estimators.adiabatic_population(output.sim, final)
    elseif output.type == :diabatic
        populations = Estimators.diabatic_population(output.sim, final)
    else
        throw(ArgumentError("$(output.type) not recognised.
            Only `:diabatic` or `:adiabatic` accepted."))
    end
    output = ComponentVector(
        reflection=zeros(NQCModels.nstates(output.sim)),
        transmission=zeros(NQCModels.nstates(output.sim))
    )
    x = DynamicsUtils.get_positions(final)[1]
    if x > 0 # If final position past 0 then we count as transmission 
        output.transmission .= populations
    else # If final position left of 0 then we count as reflection
        output.reflection .= populations
    end
    return output
end
export OutputStateResolvedScattering1D

end # module
