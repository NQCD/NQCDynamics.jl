using Unitful
using UnitfulAtomic
using LinearAlgebra: norm

abstract type AbstractOutput end

"""
$(TYPEDEF)

Output the end point of each trajectory.
"""
struct OutputFinal <: AbstractOutput end

(::OutputFinal)(sol, i) = (last(sol), false)

"""
$(TYPEDEF)

Output a 1 if the molecule has dissociated, 0 otherwise.

$(FIELDS)
"""
struct OutputDissociation{T,A} <: AbstractOutput
    "The maximum distance at which the two atoms can be considered bonded."
    distance::T
    "The indices of the two atoms in the molecule of interest."
    atom_indices::A
end
OutputDissociation(distance::Unitful.Quantity, atom_indices) = OutputDissociation(austrip(distance), atom_indices)

function (output::OutputDissociation)(sol, i)
    R = Dynamics.get_positions(last(sol))
    dissociated = norm(R[:,output.atom_indices[1]] .- R[:,output.atom_indices[2]]) > output.distance
    return dissociated ? (1, false) : (0, false)
end


"""
$(TYPEDEF)

Output the population of each diabatic state.
"""
struct OutputDiabaticPopulation <: AbstractOutput end
(output::OutputDiabaticPopulation)(sol, i) = (Dynamics.get_population.(sol.p, sol.u), false)
