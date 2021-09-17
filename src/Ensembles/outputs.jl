
using UnitfulAtomic: austrip
using LinearAlgebra: norm

using ..InitialConditions: QuantisedDiatomic
using NonadiabaticMolecularDynamics: Estimators

abstract type AbstractOutput end

"""
Output the end point of each trajectory.
"""
struct OutputFinal <: AbstractOutput end

(::OutputFinal)(sol, i) = (last(sol), false)

"""
Output a 1 if the molecule has dissociated, 0 otherwise.
"""
struct OutputDissociation{T,A} <: AbstractOutput
    "The maximum distance at which the two atoms can be considered bonded."
    distance::T
    "The indices of the two atoms in the molecule of interest."
    atom_indices::A
end
OutputDissociation(distance, atom_indices) = OutputDissociation(austrip(distance), atom_indices)

function (output::OutputDissociation)(sol, i)
    R = DynamicsMethods.get_positions(last(sol))
    dissociated = norm(R[:,output.atom_indices[1]] .- R[:,output.atom_indices[2]]) > output.distance
    return dissociated ? (1, false) : (0, false)
end


"""
Output the population of each diabatic state.
"""
struct OutputDiabaticPopulation{S} <: AbstractOutput
    sim::S
end
(output::OutputDiabaticPopulation)(sol, i) = (Estimators.diabatic_population.(output.sim, sol.u), false)

"""
Output the population of each adiabatic state.
"""
struct OutputAdiabaticPopulation{S} <: AbstractOutput
    sim::S
end
(output::OutputAdiabaticPopulation)(sol, i) = (Estimators.get_adiabatic_population.(output.sim, sol.u), false)

"""
Output the vibrational and rotational quantum numbers of the final image.
"""
struct OutputQuantisedDiatomic{S,H,V} <: AbstractOutput
    sim::S
    height::H
    normal_vector::V
end
OutputQuantisedDiatomic(sim; height=10, normal_vector=[0, 0, 1]) = OutputQuantisedDiatomic(sim, height, normal_vector)

function (output::OutputQuantisedDiatomic)(sol, i)
    final = last(sol.u) 
    ν, J = QuantisedDiatomic.quantise_diatomic(output.sim,
        DynamicsUtils.get_velocities(final), DynamicsUtils.get_positions(final);
        height=output.height, normal_vector=output.normal_vector)
    return ((ν, J), false)
end

"""
Output a matrix with `size=(n_states, 2)` with transmission/reflection probabilities for
each state.
"""
struct OutputStateResolvedScattering1D{S} <: AbstractOutput
    sim::S
end
function (output::OutputStateResolvedScattering1D)(sol, i)
    final = last(sol.u) # get final configuration from trajectory
    populations = Estimators.adiabatic_population(output.sim, final)
    output = zeros(2, 2) # Initialise output, left column reflection on state 1/2, right column transmission
    x = DynamicsUtils.get_positions(final)[1]
    if x > 0 # If final position past 0 then we count as transmission 
        output[:,2] .= populations
    else # If final position left of 0 then we count as reflection
        output[:,1] .= populations
    end
    return (output, false)
end
