module AdvancedASEFunctions
using NQCModels, PyCall


"""
This module contains methods related to the NQCModels ASE interface that need access to Python types. (e.g. constraint checking)
"""

function NQCModels.mobileatoms(model::NQCModels.AdiabaticModels.AdiabaticASEModel, n::Int)
	ase_constraints = [constraint for constraint in model.atoms.constraints]
	if length(ase_constraints) == 0
        return 1:length(model.atoms)
	else
        return symdiff(1:length(model.atoms), [constraint.get_indices() .+ 1 for constraint in ase_constraints]...)
    end
end

end
