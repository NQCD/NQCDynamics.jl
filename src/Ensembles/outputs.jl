
using UnitfulAtomic: austrip
using LinearAlgebra: norm
using ComponentArrays: ComponentVector

using ..InitialConditions: QuantisedDiatomic
using NQCDynamics: Estimators
using NQCModels: nstates
using NQCDistributions: Diabatic, Adiabatic

function output_template(output, u0)
    return zero(output(ComponentVector(u=[u0]), 1)[1])
end

"""
Output the end point of each trajectory.
"""
struct OutputFinal end

(::OutputFinal)(sol) = last(sol.u)

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

function (output::OutputDissociation)(sol, i)
    R = DynamicsUtils.get_positions(last(sol.u))
    dissociated = norm(R[:,output.atom_indices[1]] .- R[:,output.atom_indices[2]]) > output.distance
    return dissociated ? (1, false) : (0, false)
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

function (output::OutputQuantisedDiatomic)(sol, i)
    final = last(sol.u) 
    ν, J = QuantisedDiatomic.quantise_diatomic(output.sim,
        DynamicsUtils.get_velocities(final), DynamicsUtils.get_positions(final);
        height=output.height, normal_vector=output.normal_vector)
    return ((ν, J), false)
end
output_template(::OutputQuantisedDiatomic, u0) = (0.0, 0.0)

"""
Output a `ComponentVector` with fields `reflection` and `transmission` containing
the probability of the outcome.
Each index in the arrays refers to the adiabatic state.
"""
struct OutputStateResolvedScattering1D{S}
    sim::S
    type::Symbol
end
function (output::OutputStateResolvedScattering1D)(sol)
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
        reflection=zeros(nstates(output.sim)),
        transmission=zeros(nstates(output.sim))
    )
    x = DynamicsUtils.get_positions(final)[1]
    if x > 0 # If final position past 0 then we count as transmission 
        output.transmission .= populations
    else # If final position left of 0 then we count as reflection
        output.reflection .= populations
    end
    return output
end
