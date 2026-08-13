using NQCModels: ndofs

"""
    CompositeFrictionModel{M,F} <: ClassicalFrictionModel

Combine an `ClassicalModel` with an `ElectronicFrictionProvider`.

This allows for arbitrary composition of potentials and friction providers,
such that any classical model can be augmented with any form of electronic friction.
"""
struct CompositeFrictionModel{M<:ClassicalModel,F} <: ClassicalFrictionModel
    pes_model::M
    friction_model::F
end

function NQCModels.ndofs(model::CompositeFrictionModel)
    model_dofs = ndofs(model.pes_model)
    friction_dofs = ndofs(model.friction_model)

    if model_dofs == friction_dofs
        return model_dofs
    else
        throw(ArgumentError("Dimensionality of PES and friction models does not match.")) 
    end
end

function NQCModels.potential(model::CompositeFrictionModel, R::AbstractMatrix)
    return NQCModels.potential(model.pes_model, R)
end

function NQCModels.potential!(model::CompositeFrictionModel, V::AbstractMatrix, R::AbstractMatrix)
    V .= NQCModels.potential(model.pes_model, R)
end

function NQCModels.derivative!(model::CompositeFrictionModel, D, R::AbstractMatrix)
    NQCModels.derivative!(model.pes_model, D, R)
end

function friction!(model::CompositeFrictionModel, F, R)
    friction!(model.friction_model, F, R)
end
